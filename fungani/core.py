from Bio import SeqIO


def chunks(lst, n):
    if n == 0:
        n = 1
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


def blast(record, db):
    pass


def parse_results(dir):
    pass
