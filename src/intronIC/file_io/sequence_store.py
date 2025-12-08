"""
Streaming sequence storage using SQLite.

This module provides temporary storage for intron sequences during streaming mode.
Sequences are written to SQLite as they're extracted (supporting parallel workers),
then read back after scoring to write the final .introns.iic file with scores.

Design for parallel safety:
- WAL mode enables concurrent reads/writes from multiple processes
- Each worker opens its own connection (connections are not shared across processes)
- Batch inserts with single commit per contig for performance
- Main process reads all rows after all workers complete

Usage:
    # Initialize store (creates database file)
    store = StreamingSequenceStore.create(output_dir / "sequences.db")

    # In worker processes (parallel extraction)
    store = StreamingSequenceStore(db_path)
    store.insert_batch(introns_with_sequences)
    store.close()

    # In main process (after scoring)
    store = StreamingSequenceStore(db_path)
    for row in store.iter_all():
        score = score_lookup[row.intron_id]
        write_final_line(row, score)
    store.close()
    store.cleanup()  # Delete temp database

Author: intronIC refactoring project
Date: 2025-11-30
"""

import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, List, Optional, Union
from contextlib import contextmanager

from intronIC.core.intron import Intron


# Schema version for future compatibility
SCHEMA_VERSION = 1

# SQL statements
CREATE_SCHEMA = """
-- Streaming sequence storage for intronIC
-- Stores sequences temporarily during extraction, before scoring completes

-- Metadata table for schema versioning
CREATE TABLE IF NOT EXISTS metadata (
    key TEXT PRIMARY KEY,
    value TEXT
);

-- Main sequences table
-- Stores all data needed to reconstruct the .introns.iic output line
CREATE TABLE IF NOT EXISTS sequences (
    intron_id TEXT PRIMARY KEY,
    formatted_name TEXT NOT NULL,  -- Formatted output name (e.g., CaeEle-gene@transcript_N(size))
    upstream_flank TEXT,
    seq TEXT NOT NULL,
    downstream_flank TEXT,
    terminal_dnts TEXT,
    -- Preserve insertion order for deterministic output
    insertion_order INTEGER
);

-- Index not strictly needed since we read all rows, but helps with debugging
CREATE INDEX IF NOT EXISTS idx_intron_id ON sequences(intron_id);
"""

INSERT_SEQUENCE = """
INSERT INTO sequences (intron_id, formatted_name, upstream_flank, seq, downstream_flank, terminal_dnts, insertion_order)
VALUES (?, ?, ?, ?, ?, ?, ?)
"""

SELECT_ALL = """
SELECT intron_id, formatted_name, upstream_flank, seq, downstream_flank, terminal_dnts
FROM sequences
ORDER BY insertion_order
"""

SELECT_COUNT = "SELECT COUNT(*) FROM sequences"


@dataclass(frozen=True, slots=True)
class SequenceRow:
    """
    A row from the sequences table.

    This is the data structure returned when iterating over stored sequences.
    Contains all fields needed to write the final .introns.iic output line.
    """
    intron_id: str
    formatted_name: str  # Formatted output name (e.g., CaeEle-gene@transcript_N(size))
    upstream_flank: Optional[str]
    seq: str
    downstream_flank: Optional[str]
    terminal_dnts: Optional[str]


class StreamingSequenceStore:
    """
    SQLite-backed storage for sequences during streaming extraction.

    Thread/process safety:
    - Uses WAL mode for concurrent access from multiple processes
    - Each process should create its own connection (not shared)
    - Batch inserts with single commit per contig for performance

    The store is designed for a specific workflow:
    1. Main process creates the database with create()
    2. Worker processes open connections and insert batches
    3. After all workers complete, main process reads all rows
    4. Main process calls cleanup() to delete the database

    Attributes:
        db_path: Path to the SQLite database file
        conn: SQLite connection (None if not connected)
    """

    def __init__(self, db_path: Union[str, Path]):
        """
        Open a connection to an existing sequence store.

        Args:
            db_path: Path to SQLite database file (must exist)

        Note:
            Use StreamingSequenceStore.create() to create a new database.
            This constructor is for opening existing databases.
        """
        self.db_path = Path(db_path)
        self.conn: Optional[sqlite3.Connection] = None
        self._insertion_counter = 0
        self._connect()

    @classmethod
    def create(cls, db_path: Union[str, Path]) -> 'StreamingSequenceStore':
        """
        Create a new sequence store database.

        This is the factory method for creating new stores. It:
        1. Creates the database file
        2. Initializes the schema
        3. Sets up WAL mode for concurrent access
        4. Returns an open connection

        Args:
            db_path: Path where database will be created

        Returns:
            StreamingSequenceStore instance with open connection

        Raises:
            FileExistsError: If database already exists
        """
        db_path = Path(db_path)

        if db_path.exists():
            raise FileExistsError(
                f"Sequence store already exists at {db_path}. "
                "Delete it first or use a different path."
            )

        # Create parent directory if needed
        db_path.parent.mkdir(parents=True, exist_ok=True)

        # Create database and initialize schema
        conn = sqlite3.connect(str(db_path))
        try:
            # Enable WAL mode for concurrent access
            conn.execute("PRAGMA journal_mode=WAL")
            # Reduce sync for performance (we can recreate if crash)
            conn.execute("PRAGMA synchronous=NORMAL")

            # Create schema
            conn.executescript(CREATE_SCHEMA)

            # Store schema version
            conn.execute(
                "INSERT INTO metadata (key, value) VALUES (?, ?)",
                ("schema_version", str(SCHEMA_VERSION))
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
        # Enable WAL mode (idempotent)
        self.conn.execute("PRAGMA journal_mode=WAL")
        # Reduce sync for performance
        self.conn.execute("PRAGMA synchronous=NORMAL")
        # Use Row factory for named access
        self.conn.row_factory = sqlite3.Row

        # Get current max insertion_order for this connection
        cursor = self.conn.execute(
            "SELECT COALESCE(MAX(insertion_order), 0) FROM sequences"
        )
        self._insertion_counter = cursor.fetchone()[0]

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

    def insert_batch(
        self,
        introns: List[Intron],
        formatted_names: Optional[List[str]] = None
    ) -> int:
        """
        Insert a batch of introns with their sequences.

        This is the primary method for storing sequences during extraction.
        It inserts all introns in a single transaction for performance.

        Args:
            introns: List of Intron objects with sequences populated
            formatted_names: Optional list of formatted output names (parallel to introns).
                           If not provided, uses intron.intron_id as the name.

        Returns:
            Number of introns inserted

        Raises:
            ValueError: If any intron lacks sequences
            sqlite3.IntegrityError: If duplicate intron_id
        """
        if self.conn is None:
            self._connect()

        # Use intron_id as fallback for formatted name
        if formatted_names is None:
            formatted_names = [intron.intron_id for intron in introns]

        rows = []
        for intron, formatted_name in zip(introns, formatted_names):
            if intron.sequences is None or intron.sequences.seq is None:
                raise ValueError(
                    f"Intron {intron.intron_id} has no sequence. "
                    "Cannot store in sequence store."
                )

            self._insertion_counter += 1
            rows.append((
                intron.intron_id,
                formatted_name,
                intron.sequences.upstream_flank,
                intron.sequences.seq,
                intron.sequences.downstream_flank,
                intron.sequences.terminal_dinucleotides,
                self._insertion_counter
            ))

        # Batch insert with single commit
        self.conn.executemany(INSERT_SEQUENCE, rows)
        self.conn.commit()

        return len(rows)

    def insert_one(self, intron: Intron) -> None:
        """
        Insert a single intron with its sequence.

        For most use cases, prefer insert_batch() for better performance.

        Args:
            intron: Intron object with sequences populated

        Raises:
            ValueError: If intron lacks sequences
        """
        self.insert_batch([intron])

    def iter_all(self) -> Iterator[SequenceRow]:
        """
        Iterate over all stored sequences in insertion order.

        This is used after scoring to read sequences back and write
        the final .introns.iic file with scores.

        Yields:
            SequenceRow objects in insertion order

        Example:
            for row in store.iter_all():
                score = score_lookup.get(row.intron_id, "NA")
                write_line(row.intron_id, score, row.upstream_flank, ...)
        """
        if self.conn is None:
            self._connect()

        cursor = self.conn.execute(SELECT_ALL)

        for row in cursor:
            yield SequenceRow(
                intron_id=row["intron_id"],
                formatted_name=row["formatted_name"],
                upstream_flank=row["upstream_flank"],
                seq=row["seq"],
                downstream_flank=row["downstream_flank"],
                terminal_dnts=row["terminal_dnts"]
            )

    def count(self) -> int:
        """
        Get the number of stored sequences.

        Returns:
            Total count of sequences in store
        """
        if self.conn is None:
            self._connect()

        cursor = self.conn.execute(SELECT_COUNT)
        return cursor.fetchone()[0]

    def get(self, intron_id: str) -> Optional[SequenceRow]:
        """
        Get a specific sequence by intron ID.

        This is primarily for debugging/testing. For production use,
        prefer iter_all() which is more efficient for reading all rows.

        Args:
            intron_id: ID of intron to retrieve

        Returns:
            SequenceRow if found, None otherwise
        """
        if self.conn is None:
            self._connect()

        cursor = self.conn.execute(
            "SELECT intron_id, formatted_name, upstream_flank, seq, downstream_flank, terminal_dnts "
            "FROM sequences WHERE intron_id = ?",
            (intron_id,)
        )

        row = cursor.fetchone()
        if row is None:
            return None

        return SequenceRow(
            intron_id=row["intron_id"],
            formatted_name=row["formatted_name"],
            upstream_flank=row["upstream_flank"],
            seq=row["seq"],
            downstream_flank=row["downstream_flank"],
            terminal_dnts=row["terminal_dnts"]
        )

    @contextmanager
    def transaction(self):
        """
        Context manager for explicit transaction control.

        Use this when you need to insert multiple batches atomically.

        Example:
            with store.transaction():
                store.insert_batch(batch1)
                store.insert_batch(batch2)
                # Both batches committed together
        """
        if self.conn is None:
            self._connect()

        try:
            yield
            self.conn.commit()
        except Exception:
            self.conn.rollback()
            raise

    def __enter__(self) -> 'StreamingSequenceStore':
        """Support context manager protocol."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """Close connection on context exit."""
        self.close()

    def __repr__(self) -> str:
        """Developer-friendly representation."""
        status = "connected" if self.conn else "disconnected"
        return f"StreamingSequenceStore({self.db_path}, {status})"


# Module-level function for worker process initialization
_worker_store: Optional[StreamingSequenceStore] = None


def init_worker_store(db_path: str) -> None:
    """
    Initialize sequence store for a worker process.

    This is designed to be used as a Pool initializer, similar to
    init_worker_genome in indexed_genome.py.

    Args:
        db_path: Path to the SQLite database

    Example:
        with Pool(
            processes=n_processes,
            initializer=init_worker_store,
            initargs=(str(db_path),)
        ) as pool:
            pool.starmap(worker_func, inputs)
    """
    global _worker_store
    _worker_store = StreamingSequenceStore(db_path)


def get_worker_store() -> StreamingSequenceStore:
    """
    Get the sequence store for the current worker process.

    Must be called after init_worker_store().

    Returns:
        StreamingSequenceStore instance for this worker

    Raises:
        RuntimeError: If store not initialized
    """
    global _worker_store
    if _worker_store is None:
        raise RuntimeError(
            "Worker store not initialized. "
            "Call init_worker_store() first (typically via Pool initializer)."
        )
    return _worker_store


def close_worker_store() -> None:
    """
    Close the worker process sequence store.

    Call this when worker is done (e.g., in atexit handler).
    """
    global _worker_store
    if _worker_store is not None:
        _worker_store.close()
        _worker_store = None
