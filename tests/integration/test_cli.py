"""
Integration tests for CLI interface.

Tests argument parsing, configuration building, and basic CLI functionality.
"""

import shutil
import tempfile
from pathlib import Path

import pytest

from intronIC.cli.args import IntronICArgumentParser
from intronIC.cli.config import IntronICConfig, ScoringRegions
from intronIC.cli.progress import IntronICProgressReporter


class TestArgumentParser:
    """Test argument parser functionality."""

    def test_minimal_required_arguments(self, tmp_path):
        """Test parser with minimal required arguments."""
        # Create dummy files
        genome = tmp_path / "genome.fa"
        annotation = tmp_path / "annotation.gff3"
        genome.write_text(">chr1\nACTG\n")
        annotation.write_text("chr1\t.\tgene\t1\t100\t.\t+\t.\tID=gene1\n")

        parser = IntronICArgumentParser()
        args = parser.parse_args(
            ["-n", "test_species", "-g", str(genome), "-a", str(annotation)]
        )

        assert args.species_name == "test_species"
        assert args.genome == genome
        assert args.annotation == annotation

    def test_all_input_modes(self, tmp_path):
        """Test parser accepts different input modes."""
        # Create dummy files
        genome = tmp_path / "genome.fa"
        annotation = tmp_path / "annotation.gff3"
        bed = tmp_path / "introns.bed"
        sequences = tmp_path / "introns.iic"

        for f in [genome, annotation, bed, sequences]:
            f.write_text("dummy\n")

        parser = IntronICArgumentParser()

        # Annotation mode
        args1 = parser.parse_args(
            ["-n", "species", "-g", str(genome), "-a", str(annotation)]
        )
        assert args1.annotation is not None

        # BED mode
        args2 = parser.parse_args(["-n", "species", "-g", str(genome), "-b", str(bed)])
        assert args2.bed is not None

        # Sequences mode
        args3 = parser.parse_args(["-n", "species", "-q", str(sequences)])
        assert args3.sequence_file is not None

    def test_scoring_options(self, tmp_path):
        """Test scoring-related options."""
        genome = tmp_path / "genome.fa"
        annotation = tmp_path / "annotation.gff3"
        genome.write_text(">chr1\nACTG\n")
        annotation.write_text("chr1\t.\tgene\t1\t100\t.\t+\t.\tID=gene1\n")

        parser = IntronICArgumentParser()
        args = parser.parse_args(
            [
                "-n",
                "species",
                "-g",
                str(genome),
                "-a",
                str(annotation),
                "-t",
                "95",
                "--no_nc",
            ]
        )

        assert args.threshold == 95.0
        assert args.no_nc is True

    def test_performance_options(self, tmp_path):
        """Test performance-related options."""
        genome = tmp_path / "genome.fa"
        annotation = tmp_path / "annotation.gff3"
        genome.write_text(">chr1\nACTG\n")
        annotation.write_text("chr1\t.\tgene\t1\t100\t.\t+\t.\tID=gene1\n")

        parser = IntronICArgumentParser()
        args = parser.parse_args(
            [
                "-n",
                "species",
                "-g",
                str(genome),
                "-a",
                str(annotation),
                "-p",
                "8",
                "--cv_processes",
                "4",
            ]
        )

        assert args.processes == 8
        assert args.cv_processes == 4

    def test_training_options(self, tmp_path):
        """Test training-related options (train subcommand)."""
        parser = IntronICArgumentParser()
        # Training options are now under the 'train' subcommand
        args = parser.parse_args(
            [
                "train",
                "-n",
                "species",
                "-C",
                "0.1",
                "--n_models",
                "5",
                "--seed",
                "123",
            ]
        )

        assert args.command == "train"
        assert args.C == 0.1
        assert args.n_models == 5
        assert args.seed == 123

    def test_custom_scoring_regions(self, tmp_path):
        """Test custom scoring region coordinates."""
        genome = tmp_path / "genome.fa"
        annotation = tmp_path / "annotation.gff3"
        genome.write_text(">chr1\nACTG\n")
        annotation.write_text("chr1\t.\tgene\t1\t100\t.\t+\t.\tID=gene1\n")

        parser = IntronICArgumentParser()
        args = parser.parse_args(
            [
                "-n",
                "species",
                "-g",
                str(genome),
                "-a",
                str(annotation),
                "--five_score_coords",
                "-5",
                "10",
                "--bp_region_coords",
                "-60",
                "-10",
                "--three_score_coords",
                "-8",
                "5",
            ]
        )

        assert args.five_score_coords == [-5, 10]
        assert args.bp_region_coords == [-60, -10]
        assert args.three_score_coords == [-8, 5]

    def test_missing_required_argument(self):
        """Test parser rejects missing required arguments."""
        parser = IntronICArgumentParser()
        with pytest.raises(SystemExit):
            parser.parse_args([])  # Missing -n species_name

    def test_missing_input_source(self, tmp_path):
        """Test parser rejects missing input sources."""
        parser = IntronICArgumentParser()
        with pytest.raises(SystemExit):
            parser.parse_args(["-n", "species"])  # No input files

    def test_invalid_threshold(self, tmp_path):
        """Test parser rejects invalid threshold values."""
        genome = tmp_path / "genome.fa"
        annotation = tmp_path / "annotation.gff3"
        genome.write_text(">chr1\nACTG\n")
        annotation.write_text("chr1\t.\tgene\t1\t100\t.\t+\t.\tID=gene1\n")

        parser = IntronICArgumentParser()

        with pytest.raises(SystemExit):
            parser.parse_args(
                [
                    "-n",
                    "species",
                    "-g",
                    str(genome),
                    "-a",
                    str(annotation),
                    "-t",
                    "150",  # Invalid: > 100
                ]
            )

    def test_nonexistent_file(self, tmp_path):
        """Test parser rejects nonexistent files."""
        parser = IntronICArgumentParser()
        with pytest.raises(SystemExit):
            parser.parse_args(
                ["-n", "species", "-g", "nonexistent.fa", "-a", "nonexistent.gff3"]
            )


class TestConfiguration:
    """Test configuration building from arguments."""

    def test_config_from_annotation_args(self, tmp_path):
        """Test building config from annotation mode arguments."""
        genome = tmp_path / "genome.fa"
        annotation = tmp_path / "annotation.gff3"
        genome.write_text(">chr1\nACTG\n")
        annotation.write_text("chr1\t.\tgene\t1\t100\t.\t+\t.\tID=gene1\n")

        parser = IntronICArgumentParser()
        args = parser.parse_args(
            [
                "-n",
                "test_species",
                "-g",
                str(genome),
                "-a",
                str(annotation),
                "-t",
                "95",
                "-p",
                "4",
            ]
        )

        config = IntronICConfig.from_args(args)

        assert config.input.genome == genome
        assert config.input.annotation == annotation
        assert config.input.mode == "annotation"
        assert config.scoring.threshold == 95.0
        assert config.performance.processes == 4
        assert config.output.species_name == "test_species"

    def test_config_from_bed_args(self, tmp_path):
        """Test building config from BED mode arguments."""
        genome = tmp_path / "genome.fa"
        bed = tmp_path / "introns.bed"
        genome.write_text(">chr1\nACTG\n")
        bed.write_text("chr1\t100\t200\tintron1\t0\t+\n")

        parser = IntronICArgumentParser()
        args = parser.parse_args(
            ["-n", "test_species", "-g", str(genome), "-b", str(bed)]
        )

        config = IntronICConfig.from_args(args)

        assert config.input.genome == genome
        assert config.input.bed == bed
        assert config.input.mode == "bed"

    def test_config_from_sequences_args(self, tmp_path):
        """Test building config from sequences mode arguments."""
        sequences = tmp_path / "introns.iic"
        sequences.write_text("intron1\tACTG\tGTAAG\tTTCAG\n")

        parser = IntronICArgumentParser()
        args = parser.parse_args(["-n", "test_species", "-q", str(sequences)])

        config = IntronICConfig.from_args(args)

        assert config.input.sequence_file == sequences
        assert config.input.mode == "sequences"

    def test_config_scoring_regions(self, tmp_path):
        """Test configuration with custom scoring regions."""
        genome = tmp_path / "genome.fa"
        annotation = tmp_path / "annotation.gff3"
        genome.write_text(">chr1\nACTG\n")
        annotation.write_text("chr1\t.\tgene\t1\t100\t.\t+\t.\tID=gene1\n")

        parser = IntronICArgumentParser()
        args = parser.parse_args(
            [
                "-n",
                "species",
                "-g",
                str(genome),
                "-a",
                str(annotation),
                "--five_score_coords",
                "-5",
                "10",
                "--bp_region_coords",
                "-60",
                "-10",
                "--three_score_coords",
                "-8",
                "5",
            ]
        )

        config = IntronICConfig.from_args(args)

        assert config.scoring.scoring_regions.five_start == -5
        assert config.scoring.scoring_regions.five_end == 10
        assert config.scoring.scoring_regions.bp_start == -60
        assert config.scoring.scoring_regions.bp_end == -10
        assert config.scoring.scoring_regions.three_start == -8
        assert config.scoring.scoring_regions.three_end == 5

    def test_config_extraction_options(self, tmp_path):
        """Test configuration with extraction options."""
        genome = tmp_path / "genome.fa"
        annotation = tmp_path / "annotation.gff3"
        genome.write_text(">chr1\nACTG\n")
        annotation.write_text("chr1\t.\tgene\t1\t100\t.\t+\t.\tID=gene1\n")

        parser = IntronICArgumentParser()
        args = parser.parse_args(
            [
                "-n",
                "species",
                "-g",
                str(genome),
                "-a",
                str(annotation),
                "-f",
                "cds",
                "--min_intron_len",
                "50",
                "--flank_len",
                "100",
                "-i",
                "-v",
            ]
        )

        config = IntronICConfig.from_args(args)

        assert config.extraction.feature_type == "cds"
        assert config.extraction.min_intron_len == 50
        assert config.extraction.flank_len == 100
        assert config.extraction.allow_multiple_isoforms is True
        assert config.extraction.no_intron_overlap is True

    def test_config_training_options(self, tmp_path):
        """Test configuration with training options (train subcommand)."""
        parser = IntronICArgumentParser()
        # Training options are now under the 'train' subcommand
        args = parser.parse_args(
            [
                "train",
                "-n",
                "species",
                "-C",
                "0.5",
                "--n_models",
                "3",
                "--seed",
                "999",
            ]
        )

        config = IntronICConfig.from_args(args)

        assert config.training.fixed_C == 0.5
        assert config.training.n_models == 3
        assert config.training.seed == 999

    def test_config_output_paths(self, tmp_path):
        """Test configuration generates correct output paths."""
        genome = tmp_path / "genome.fa"
        annotation = tmp_path / "annotation.gff3"
        genome.write_text(">chr1\nACTG\n")
        annotation.write_text("chr1\t.\tgene\t1\t100\t.\t+\t.\tID=gene1\n")

        parser = IntronICArgumentParser()
        args = parser.parse_args(
            [
                "-n",
                "homo_sapiens",
                "-g",
                str(genome),
                "-a",
                str(annotation),
                "-o",
                str(tmp_path),
            ]
        )

        config = IntronICConfig.from_args(args)

        assert config.output.base_filename == "homo_sapiens"
        assert (
            config.output.get_output_path(".meta.iic")
            == tmp_path / "homo_sapiens.meta.iic"
        )
        assert (
            config.output.get_output_path(".bed.iic")
            == tmp_path / "homo_sapiens.bed.iic"
        )


class TestProgressReporter:
    """Test progress reporter functionality."""

    def test_reporter_initialization(self):
        """Test reporter can be initialized."""
        reporter = IntronICProgressReporter(quiet=False)
        assert reporter.console is not None
        assert reporter.quiet is False

    def test_quiet_mode(self):
        """Test quiet mode suppresses output."""
        reporter = IntronICProgressReporter(quiet=True)
        assert reporter.quiet is True

        # These should not raise errors even in quiet mode
        reporter.print_info("Test message")
        reporter.print_success("Success message")
        reporter.print_warning("Warning message")

    def test_classification_summary(self):
        """Test classification summary rendering."""
        reporter = IntronICProgressReporter(quiet=False)

        # Should not raise errors
        reporter.print_classification_summary(
            total=1000, u12_count=5, u2_count=995, threshold=90.0
        )

    def test_file_tree(self):
        """Test file tree rendering."""
        reporter = IntronICProgressReporter(quiet=False)

        output_files = {
            "Metadata": "/path/to/output.meta.iic",
            "BED": "/path/to/output.bed.iic",
            "Sequences": "/path/to/output.seqs.iic",
        }

        # Should not raise errors
        reporter.print_file_tree(output_files)

    def test_pipeline_steps(self):
        """Test pipeline steps rendering."""
        reporter = IntronICProgressReporter(quiet=False)

        steps = ["Load input data", "Extract introns", "Score introns", "Classify"]

        # Should not raise errors
        reporter.print_pipeline_steps(steps)
        reporter.print_pipeline_steps(steps, current_step=2)

    def test_stats_table(self):
        """Test statistics table rendering."""
        reporter = IntronICProgressReporter(quiet=False)

        stats = {
            "Total introns": 10000,
            "U12-type": 50,
            "U2-type": 9950,
            "U12 percentage": 0.005,
            "Average SVM score": 85.3,
        }

        # Should not raise errors
        reporter.print_stats_table(stats, title="Pipeline Statistics")


class TestTrueStreamingClassification:
    """Tests for true streaming per-contig classification mode."""

    def test_true_streaming_requires_pretrained_model(self):
        """Test that true streaming mode requires a pre-trained model."""
        from unittest.mock import MagicMock, PropertyMock

        from intronIC.cli.config import InputConfig, IntronICConfig, TrainingConfig
        from intronIC.cli.main import classify_streaming_per_contig
        from intronIC.cli.messenger import UnifiedMessenger
        from intronIC.cli.progress import IntronICProgressReporter

        # Create config without pretrained model using proper nested mocks
        config = MagicMock()
        config.training = MagicMock(spec=TrainingConfig)
        config.training.pretrained_model_path = None
        config.input = MagicMock(spec=InputConfig)
        config.input.mode = "annotation"

        messenger = MagicMock(spec=UnifiedMessenger)
        reporter = IntronICProgressReporter(quiet=True)

        with pytest.raises(ValueError, match="pre-trained model"):
            classify_streaming_per_contig(config, messenger, reporter)

    def test_true_streaming_requires_annotation_mode(self):
        """Test that true streaming mode only works with annotation input."""
        from pathlib import Path
        from unittest.mock import MagicMock

        from intronIC.cli.config import InputConfig, IntronICConfig, TrainingConfig
        from intronIC.cli.main import classify_streaming_per_contig
        from intronIC.cli.messenger import UnifiedMessenger
        from intronIC.cli.progress import IntronICProgressReporter

        # Create config with pretrained model but BED mode using proper nested mocks
        config = MagicMock()
        config.training = MagicMock(spec=TrainingConfig)
        config.training.pretrained_model_path = Path("/fake/model.pkl")
        config.input = MagicMock(spec=InputConfig)
        config.input.mode = "bed"

        messenger = MagicMock(spec=UnifiedMessenger)
        reporter = IntronICProgressReporter(quiet=True)

        with pytest.raises(ValueError, match="annotation input mode"):
            classify_streaming_per_contig(config, messenger, reporter)
