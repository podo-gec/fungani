from fungani.cli import parse_args


def test_parser_default_outdir():
    parser = parse_args(["reference.fasta", "test.fasta"])
    assert parser.outdir is None


def test_parser_user_outdir():
    parser = parse_args(["reference.fasta", "test.fasta", "-o", "tmp"])
    assert parser.outdir is not None
