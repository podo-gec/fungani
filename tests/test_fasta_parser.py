import os
import tempfile

from fqfa.fasta import fasta


DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), "data")
A = os.path.join(DATA_DIR, "aa.fasta")
B = os.path.join(DATA_DIR, "bb.fasta")
C = os.path.join(DATA_DIR, "random_its.fasta")


def test_fasta_reader():
    with open(C, "r") as file:
        records = list(fasta.parse_fasta_records(file))
    assert len(records) == 1683


def test_fasta_writer():
    # NOTE: This may fail on Windows
    with open(C, "r") as file:
        records = list(fasta.parse_fasta_records(file))

    tmpfile = tempfile.NamedTemporaryFile(delete=False)
    try:
        with open(tmpfile.name, "w") as outfile:
            for record in records:
                fasta.write_fasta_record(outfile, record[0], record[1])
        with open(tmpfile.name, "r") as infile:
            data = list(fasta.parse_fasta_records(infile))
        assert len(data) == 1683
    finally:
        tmpfile.close()
        os.unlink(tmpfile.name)
