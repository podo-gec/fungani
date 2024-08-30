import os


LOG_FILE = os.path.join(os.path.expanduser("~"), "fungani.log")
OUTFILE_FWD = os.path.join(os.path.expanduser("~"), "fungani_fwd.csv")
OUTFILE_REV = os.path.join(os.path.expanduser("~"), "fungani_rev.csv")


def test_log_file_exists():
    assert os.path.isfile(LOG_FILE)


def test_result_file_exists():
    assert os.path.isfile(OUTFILE_FWD)
    assert os.path.isfile(OUTFILE_REV)
