import pytest
import os
import glob
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from src.ab1_helpers import (
    ab1_reader,
    mult_ab1_to_dict,
    mult_ab1_to_single_fasta,
    blast_seq,
    blast_seqdict,
    xml_to_table,
)

SAMPLE_AB1 = "seq_data/old_ab1s/2009-05-13_A10_AAPAA406T7.ab1"
AB1_SUBFOLDER = "old_ab1s"
EXISTING_OUTPUT = "python_outputs/old_ab1s"


# --- ab1_reader ---

class TestAb1Reader:
    def test_returns_seqrecord(self):
        record = ab1_reader(SAMPLE_AB1)
        assert isinstance(record, SeqRecord)

    def test_sequence_is_nonempty(self):
        record = ab1_reader(SAMPLE_AB1)
        assert len(record.seq) > 0


# --- mult_ab1_to_dict ---

class TestMultAb1ToDict:
    def test_returns_dict(self):
        result = mult_ab1_to_dict(AB1_SUBFOLDER)
        assert isinstance(result, dict)

    def test_correct_count(self):
        result = mult_ab1_to_dict(AB1_SUBFOLDER)
        assert len(result) == 10

    def test_values_are_seq_objects(self):
        result = mult_ab1_to_dict(AB1_SUBFOLDER)
        for val in result.values():
            assert isinstance(val, Seq)

    def test_nonexistent_subfolder_returns_none(self):
        result = mult_ab1_to_dict("does_not_exist")
        assert result is None or len(result) == 0


# --- mult_ab1_to_single_fasta ---

class TestMultAb1ToSingleFasta:
    @pytest.fixture(autouse=True)
    def setup_tmp_ab1(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        os.makedirs("seq_data/old_ab1s")
        os.makedirs("python_outputs")
        src = os.path.join(os.path.dirname(__file__), "..", SAMPLE_AB1)
        import shutil
        shutil.copy(os.path.abspath(src), "seq_data/old_ab1s/")

    def test_creates_fasta_file(self):
        seqdict, output_folder = mult_ab1_to_single_fasta("old_ab1s", "test_output.fasta")
        fasta_path = os.path.join(output_folder, "test_output.fasta")

        assert os.path.exists(fasta_path)
        assert isinstance(seqdict, dict)
        assert len(seqdict) == 1

    def test_adds_fasta_extension(self):
        seqdict, output_folder = mult_ab1_to_single_fasta("old_ab1s", "no_extension")
        assert os.path.exists(os.path.join(output_folder, "no_extension.fasta"))

    def test_fasta_has_correct_header_count(self):
        seqdict, output_folder = mult_ab1_to_single_fasta("old_ab1s")
        fasta_path = os.path.join(output_folder, "my_ab1s.fasta")
        with open(fasta_path) as f:
            headers = [line for line in f if line.startswith(">")]
        assert len(headers) == len(seqdict)


# --- xml_to_table ---

class TestXmlToTable:
    def test_creates_csv(self, tmp_path):
        seqdict = mult_ab1_to_dict(AB1_SUBFOLDER)
        # Copy existing XMLs to tmp dir so we don't overwrite the real CSV
        import shutil
        tmp_output = str(tmp_path / "output")
        shutil.copytree(EXISTING_OUTPUT, tmp_output)

        xml_to_table(seqdict, tmp_output)
        csvfile = os.path.join(tmp_output, "summary_of_blast.csv")
        assert os.path.exists(csvfile)

    def test_csv_has_correct_header(self, tmp_path):
        seqdict = mult_ab1_to_dict(AB1_SUBFOLDER)
        import shutil
        tmp_output = str(tmp_path / "output")
        shutil.copytree(EXISTING_OUTPUT, tmp_output)

        xml_to_table(seqdict, tmp_output)
        csvfile = os.path.join(tmp_output, "summary_of_blast.csv")
        with open(csvfile) as f:
            header = f.readline().strip()
        assert "Sample_ID" in header
        assert "Identity(%)" in header
        assert "NCBI_Link" in header

    def test_csv_row_count_matches_xmls(self, tmp_path):
        seqdict = mult_ab1_to_dict(AB1_SUBFOLDER)
        import shutil
        tmp_output = str(tmp_path / "output")
        shutil.copytree(EXISTING_OUTPUT, tmp_output)

        xml_to_table(seqdict, tmp_output)
        csvfile = os.path.join(tmp_output, "summary_of_blast.csv")
        xml_count = len(glob.glob(os.path.join(tmp_output, "*.xml")))
        with open(csvfile) as f:
            row_count = sum(1 for _ in f) - 1  # subtract header
        assert row_count == xml_count


# --- NCBI network tests (skip in CI) ---

@pytest.mark.slow
class TestBlastSeq:
    def test_returns_readable_result(self):
        seqdict = mult_ab1_to_dict(AB1_SUBFOLDER)
        first_seq = next(iter(seqdict.values()))
        result = blast_seq(first_seq)
        content = result.read()
        assert "BlastOutput" in content


@pytest.mark.slow
class TestBlastSeqdict:
    def test_saves_xml_files(self, tmp_path):
        seqdict = mult_ab1_to_dict(AB1_SUBFOLDER)
        # Use only 1 sequence to keep it fast
        single = {k: v for k, v in list(seqdict.items())[:1]}
        output = str(tmp_path)
        blast_seqdict(single, output)
        xmls = glob.glob(os.path.join(output, "*.xml"))
        assert len(xmls) == 1
