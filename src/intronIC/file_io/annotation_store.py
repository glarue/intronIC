"""
Streaming annotation storage using SQLite.

This module provides temporary storage for parsed annotations during streaming mode.
Annotations are parsed once and stored in SQLite indexed by contig, enabling
memory-efficient per-contig retrieval without re-scanning the source file.

Design goals:
- Single parse pass: Parse annotation file once, store parsed objects
- O(1) contig access: Retrieve annotations for any contig via index
- Minimal RAM: Annotations stored on disk, only current contig in memory
- Clean interface: Hide SQLite complexity from callers

Usage:
    # Create store from annotation file (parses once)
    store = StreamingAnnotationStore.create_from_file(
        annotation_path="annotation.gff3",
        db_path=Path("temp/annotations.db"),
        parser=BioGLAnnotationParser(),
    )

    # Get list of contigs
    contigs = store.get_contigs()

    # Process each contig
    for contig in contigs:
        annotations = store.get_annotations_for_contig(contig)
        # ... process annotations ...

    # Cleanup when done
    store.cleanup()

Author: intronIC refactoring project
Date: 2025-12-02
"""

import json
import sqlite3
from pathlib import Path
from typing import Iterator, List, Optional, Union

from intronIC.file_io.parsers import AnnotationLine, BioGLAnnotationParser

# Schema version for future compatibility
ANNOTATION_SCHEMA_VERSION = 1

# SQL statements
CREATE_ANNOTATION_SCHEMA = """
-- Streaming annotation storage for intronIC
-- Stores parsed annotations indexed by contig for efficient per-contig retrieval

-- Metadata table for schema versioning
CREATE TABLE IF NOT EXISTS metadata (
    key TEXT PRIMARY KEY,
    value TEXT
);

-- Main annotations table
-- Stores all fields from AnnotationLine dataclass
CREATE TABLE IF NOT EXISTS annotations (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    contig TEXT NOT NULL,
    name TEXT NOT NULL,
    feat_type TEXT NOT NULL,
    parent_json TEXT,  -- JSON array of parent IDs
    grandparent TEXT,
    strand TEXT NOT NULL,
    start INTEGER NOT NULL,
    stop INTEGER NOT NULL,
    line_number INTEGER NOT NULL,
    phase INTEGER,
    attributes_json TEXT  -- JSON dict for additional attributes
);

-- Index for efficient per-contig retrieval
CREATE INDEX IF NOT EXISTS idx_annotations_contig ON annotations(contig);

-- Table to track unique contigs (for fast listing)
CREATE TABLE IF NOT EXISTS contigs (
    name TEXT PRIMARY KEY,
    annotation_count INTEGER DEFAULT 0
);
"""

INSERT_ANNOTATION = """
INSERT INTO annotations (
    contig, name, feat_type, parent_json, grandparent, 
    strand, start, stop, line_number, phase, attributes_json
)
VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
"""

INSERT_OR_UPDATE_CONTIG = """
INSERT INTO contigs (name, annotation_count) VALUES (?, 1)
ON CONFLICT(name) DO UPDATE SET annotation_count = annotation_count + 1
"""

SELECT_CONTIGS = "SELECT name FROM contigs ORDER BY name"

SELECT_CONTIGS_WITH_COUNTS = "SELECT name, annotation_count FROM contigs ORDER BY name"

SELECT_ANNOTATIONS_FOR_CONTIG = """
SELECT name, feat_type, parent_json, grandparent, contig, strand, 
       start, stop, line_number, phase, attributes_json
FROM annotations
WHERE contig = ?
ORDER BY line_number
"""

SELECT_CONTIG_COUNT = "SELECT COUNT(*) FROM contigs"

SELECT_TOTAL_ANNOTATIONS = "SELECT COUNT(*) FROM annotations"


class StreamingAnnotationStore:
    """
    SQLite-backed storage for annotations during streaming mode.

    This class enables memory-efficient processing of large annotation files
    by storing parsed annotations in SQLite indexed by contig. This allows:

    1. Single-pass parsing: Parse the annotation file only once
    2. Per-contig retrieval: Load only current contig's annotations into RAM
    3. Arbitrary contig counts: Handles genomes with thousands of scaffolds

    The store is designed for a specific workflow:
    1. Create store from annotation file with create_from_file()
    2. Get list of contigs with get_contigs()
    3. For each contig, get_annotations_for_contig() returns AnnotationLine objects
    4. Call cleanup() when done to delete the database

    Attributes:
        db_path: Path to the SQLite database file
        conn: SQLite connection (None if not connected)
    """

    def __init__(self, db_path: Union[str, Path]):
        """
        Open a connection to an existing annotation store.

        Args:
            db_path: Path to SQLite database file (must exist)

        Note:
            Use StreamingAnnotationStore.create_from_file() to create a new store.
            This constructor is for opening existing databases.
        """
        self.db_path = Path(db_path)
        self.conn: Optional[sqlite3.Connection] = None
        self._connect()

    @classmethod
    def create_from_file(
        cls,
        annotation_path: Union[str, Path],
        db_path: Union[str, Path],
        parser: Optional[BioGLAnnotationParser] = None,
        batch_size: int = 10000,
    ) -> "StreamingAnnotationStore":
        """
        Create a new annotation store by parsing an annotation file.

        This is the primary factory method. It:
        1. Creates the SQLite database
        2. Parses the entire annotation file
        3. Inserts all annotations into the database
        4. Returns an open connection for querying

        Args:
            annotation_path: Path to GFF3/GTF annotation file
            db_path: Path where SQLite database will be created
            parser: Optional parser instance (creates default if not provided)
            batch_size: Number of annotations to batch before committing

        Returns:
            StreamingAnnotationStore instance ready for querying

        Raises:
            FileExistsError: If database already exists
            FileNotFoundError: If annotation file doesn't exist
        """
        annotation_path = Path(annotation_path)
        db_path = Path(db_path)

        if not annotation_path.exists():
            raise FileNotFoundError(f"Annotation file not found: {annotation_path}")

        if db_path.exists():
            raise FileExistsError(
                f"Annotation store already exists at {db_path}. "
                "Delete it first or use a different path."
            )

        # Create parent directory if needed
        db_path.parent.mkdir(parents=True, exist_ok=True)

        # Use default parser if not provided
        if parser is None:
            parser = BioGLAnnotationParser(clean_names=True)

        # Create database and schema
        conn = sqlite3.connect(str(db_path))
        try:
            # Performance settings for bulk insert
            conn.execute("PRAGMA journal_mode=WAL")
            conn.execute("PRAGMA synchronous=NORMAL")
            conn.execute("PRAGMA cache_size=-64000")  # 64MB cache

            # Create schema
            conn.executescript(CREATE_ANNOTATION_SCHEMA)

            # Store schema version
            conn.execute(
                "INSERT INTO metadata (key, value) VALUES (?, ?)",
                ("schema_version", str(ANNOTATION_SCHEMA_VERSION)),
            )

            # Parse and insert annotations in batches
            batch = []
            contig_counts: dict[str, int] = {}

            for ann_line in parser.parse_file(str(annotation_path)):
                # Serialize complex fields to JSON
                parent_json = json.dumps(ann_line.parent) if ann_line.parent else None
                attributes_json = (
                    json.dumps(ann_line.attributes) if ann_line.attributes else None
                )

                batch.append(
                    (
                        ann_line.region,  # contig
                        ann_line.name,
                        ann_line.feat_type,
                        parent_json,
                        ann_line.grandparent,
                        ann_line.strand,
                        ann_line.start,
                        ann_line.stop,
                        ann_line.line_number,
                        ann_line.phase,
                        attributes_json,
                    )
                )

                # Track contig counts
                contig_counts[ann_line.region] = (
                    contig_counts.get(ann_line.region, 0) + 1
                )

                # Batch commit for performance
                if len(batch) >= batch_size:
                    conn.executemany(INSERT_ANNOTATION, batch)
                    batch = []

            # Insert remaining batch
            if batch:
                conn.executemany(INSERT_ANNOTATION, batch)

            # Insert contig summary
            conn.executemany(
                "INSERT INTO contigs (name, annotation_count) VALUES (?, ?)",
                list(contig_counts.items()),
            )

            conn.commit()

        finally:
            conn.close()

        # Return new instance with open connection
        return cls(db_path)

    def _connect(self) -> None:
        """Open database connection with appropriate settings."""
        if self.conn is not None:
            return

        self.conn = sqlite3.connect(str(self.db_path))
        # WAL mode for better read performance
        self.conn.execute("PRAGMA journal_mode=WAL")
        # Use Row factory for named access
        self.conn.row_factory = sqlite3.Row

    def close(self) -> None:
        """Close database connection."""
        if self.conn is not None:
            self.conn.close()
            self.conn = None

    def cleanup(self) -> None:
        """
        Delete the database file and associated WAL files.

        Call this after all processing is complete to clean up
        temporary files.
        """
        self.close()

        # Delete main database file
        if self.db_path.exists():
            self.db_path.unlink()

        # Delete WAL and SHM files if they exist
        wal_path = self.db_path.with_suffix(self.db_path.suffix + "-wal")
        shm_path = self.db_path.with_suffix(self.db_path.suffix + "-shm")

        if wal_path.exists():
            wal_path.unlink()
        if shm_path.exists():
            shm_path.unlink()

    def get_contigs(self) -> List[str]:
        """
        Get list of all contigs in the annotation file.

        Returns:
            Sorted list of contig names
        """
        if self.conn is None:
            self._connect()
        assert self.conn is not None  # For type checker

        cursor = self.conn.execute(SELECT_CONTIGS)
        return [row[0] for row in cursor.fetchall()]

    def get_contigs_with_counts(self) -> List[tuple[str, int]]:
        """
        Get list of all contigs with their annotation counts.

        Returns:
            List of (contig_name, annotation_count) tuples, sorted by name
        """
        if self.conn is None:
            self._connect()
        assert self.conn is not None  # For type checker

        cursor = self.conn.execute(SELECT_CONTIGS_WITH_COUNTS)
        return [(row[0], row[1]) for row in cursor.fetchall()]

    def get_contig_count(self) -> int:
        """
        Get the number of unique contigs.

        Returns:
            Number of contigs
        """
        if self.conn is None:
            self._connect()
        assert self.conn is not None  # For type checker

        cursor = self.conn.execute(SELECT_CONTIG_COUNT)
        return cursor.fetchone()[0]

    def get_total_annotations(self) -> int:
        """
        Get the total number of annotations across all contigs.

        Returns:
            Total annotation count
        """
        if self.conn is None:
            self._connect()
        assert self.conn is not None  # For type checker

        cursor = self.conn.execute(SELECT_TOTAL_ANNOTATIONS)
        return cursor.fetchone()[0]

    def get_annotations_for_contig(self, contig: str) -> List[AnnotationLine]:
        """
        Retrieve all annotations for a specific contig.

        Args:
            contig: Contig/chromosome name

        Returns:
            List of AnnotationLine objects for this contig,
            ordered by line_number (original file order)
        """
        if self.conn is None:
            self._connect()
        assert self.conn is not None  # For type checker

        cursor = self.conn.execute(SELECT_ANNOTATIONS_FOR_CONTIG, (contig,))

        annotations = []
        for row in cursor:
            # Deserialize JSON fields
            parent = json.loads(row["parent_json"]) if row["parent_json"] else []
            attributes = (
                json.loads(row["attributes_json"]) if row["attributes_json"] else {}
            )

            ann = AnnotationLine(
                name=row["name"],
                feat_type=row["feat_type"],
                parent=parent,
                grandparent=row["grandparent"],
                region=row["contig"],
                strand=row["strand"],
                start=row["start"],
                stop=row["stop"],
                line_number=row["line_number"],
                phase=row["phase"],
                attributes=attributes,
            )
            annotations.append(ann)

        return annotations

    def iter_annotations_for_contig(self, contig: str) -> Iterator[AnnotationLine]:
        """
        Iterate over annotations for a specific contig.

        Memory-efficient alternative to get_annotations_for_contig()
        for very large contigs.

        Args:
            contig: Contig/chromosome name

        Yields:
            AnnotationLine objects for this contig
        """
        if self.conn is None:
            self._connect()
        assert self.conn is not None  # For type checker

        cursor = self.conn.execute(SELECT_ANNOTATIONS_FOR_CONTIG, (contig,))

        for row in cursor:
            parent = json.loads(row["parent_json"]) if row["parent_json"] else []
            attributes = (
                json.loads(row["attributes_json"]) if row["attributes_json"] else {}
            )

            yield AnnotationLine(
                name=row["name"],
                feat_type=row["feat_type"],
                parent=parent,
                grandparent=row["grandparent"],
                region=row["contig"],
                strand=row["strand"],
                start=row["start"],
                stop=row["stop"],
                line_number=row["line_number"],
                phase=row["phase"],
                attributes=attributes,
            )

    def __enter__(self) -> "StreamingAnnotationStore":
        """Support context manager protocol."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """Close connection on context exit."""
        self.close()

    def __repr__(self) -> str:
        """Developer-friendly representation."""
        status = "connected" if self.conn else "disconnected"
        return f"StreamingAnnotationStore({self.db_path}, {status})"
