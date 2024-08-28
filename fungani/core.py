import logging
import mmap
import multiprocessing
import os
import random
import re
import shutil
import subprocess
import tempfile

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)
BLASTN = shutil.which("blastn")
MAKEBLASTDB = shutil.which("makeblastdb")


# XXX no longer necessary
def write_sequence_length(filename, pathname):
    """Write tab-delimited sequence names and sizes."""
    records = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))

    outfile = f"{filename}.sizes"
    with open(os.path.join(pathname, outfile), "w") as file:
        for key in records.items():
            file.write("\t".join([key[0], str(len(key[1].seq))]) + "\n")

    return outfile


def make_windows(pathname, args):
    splices = []
    idx = 0
    outfile = os.path.join(pathname, f"{args.test}_split.fas")
    records = SeqIO.to_dict(SeqIO.parse(args.test, "fasta"))

    for key in records.items():
        idx += 1
        end = len(key[1].seq) - args.size + 1
        for start in range(0, end, args.overlap):
            id = ":".join((str(idx), "-".join((str(start), str(start + args.size)))))
            seq = key[1].seq[start : start + args.size]
            splices.append(SeqRecord(seq, id=id, description=""))

    SeqIO.write(splices, outfile, "fasta")

    return outfile


def chunks(lst, n):
    if n == 0:
        n = 1
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


def blast(record, blast_dir, query_dir, db):
    for r in record:
        key = str(r.id)
        query = query_dir + key + ".fasta"
        SeqIO.write(r, query, "fasta")
        cmd = (
            str(BLASTN),
            "-out",
            os.path.join(blast_dir, key),
            "-outfmt",
            "0",
            "-db",
            db,
            "-query",
            query,
            "-max_target_seqs",
            "1",
        )
        subprocess.run(cmd)


def parse_results(dir):
    """Read % identity from Blast result files."""
    out, index = [], []

    for f in sorted(os.listdir(path=dir)):
        index.append(f)
        curr = os.path.join(dir, f)
        with open(curr, "r+") as ff:
            if os.stat(curr).st_size == 0:
                print(curr, "is empty")
            if os.stat(curr).st_size > 0:
                filemap = mmap.mmap(ff.fileno(), 0)
                query = re.search(rb"Identities = (\d+/\d+)", filemap)
                if query:
                    value = query.group(1).decode("utf-8")
                    rc = list(map(int, value.split("/")))
                    out.append(rc[0] / rc[1])
                else:
                    out.append(0)

    return out, index


def make_blast_db(filename):
    """Generate Blast database."""
    outfile = os.path.splitext(filename)[0] + "-db"
    cmd = (str(MAKEBLASTDB), "-dbtype", "nucl", "-in", filename, "-out", outfile)
    subprocess.run(cmd)

    return outfile


def main(args):
    logger.info("creating temporary directories")
    if args.outdir is None:
        temp_dir = tempfile.TemporaryDirectory()
        fs_tmp = temp_dir.name
    else:
        fs_tmp = args.outdir
    fs_ani_queries = os.path.join(fs_tmp, "ani_q")
    os.makedirs(fs_ani_queries)
    fs_ani_blast = os.path.join(fs_tmp, "ani_tmp")
    os.makedirs(fs_ani_blast)

    logger.info("writing Blast db")
    reference_db = make_blast_db(args.reference)

    logger.info("preparing sliced genome")
    fsplit = make_windows(fs_tmp, args)

    records = list(SeqIO.parse(fsplit, "fasta"))
    nsample = int(len(records) * args.percent / 100)
    records = random.sample(records, nsample)
    nc = int(len(records) / args.cpus)
    data = chunks(records, nc)
    p = multiprocessing.Pool(processes=args.cpus)

    logger.info("running Blast")
    result = [
        p.apply_async(
            blast,
            args=(list(x),),
            kwds={"query_dir": fs_ani_queries, "db": reference_db},
        )
        for x in data
    ]
    result = [item.get() for item in result]

    logger.info("parsing Blast")
    out, index = parse_results(fs_ani_blast)
    outfile = os.path.join(os.path.expanduser("~"), "ani_identities.txt")
    file = open(outfile, "w")
    for o in out:
        file.write(str(o) + "\n")
    file.close()
    outfile = os.path.join(os.path.expanduser("~"), "ani_coordinates.txt")
    file = open(outfile, "w")
    for i in index:
        file.write(str(i) + "\n")
    file.close()

    logger.info("deleting temporary directories")
    os.rmdir(fs_tmp)
