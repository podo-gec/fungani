import os
import fastaparser
import tempfile


DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), "data")
A = os.path.join(DATA_DIR, "aa.fasta")
B = os.path.join(DATA_DIR, "bb.fasta")
C = os.path.join(DATA_DIR, "random_its.fasta")


def test_fasta_reader():
    with open(C, "r") as file:
        data = list(fastaparser.Reader(file, parse_method="quick"))
    assert len(data) == 1683


def test_fasta_writer():
    # NOTE: This may fail on Windows
    with open(C, "r") as file:
        data = list(fastaparser.Reader(file, parse_method="quick"))

    tmpfile = tempfile.NamedTemporaryFile(delete=False)
    try:
        with open(tmpfile.name, "w") as outfile:
            writer = fastaparser.Writer(outfile)
            for d in data:
                writer.writefasta(d)
        with open(tmpfile.name, "r") as infile:
            data = list(fastaparser.Reader(infile, parse_method="quick"))
        assert len(data) == 1683
    finally:
        tmpfile.close()
        os.unlink(tmpfile.name)
