import logging
import mmap
import multiprocessing
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
from tqdm import tqdm

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)
BLASTN = shutil.which("blastn")
MAKEBLASTDB = shutil.which("makeblastdb")


def run_async(func, arglist, kwds, cpus):
    pool = multiprocessing.Pool(processes=cpus)
    jobs = [
        pool.apply_async(
            func=func,
            args=(list(arg),),
            kwds=kwds,
        )
        for arg in arglist
    ]
    pool.close()
    result = []
    for job in tqdm(
        jobs,
        total=float("inf"),
        desc=f"Computing {len(jobs)} batch (blastn) on {cpus} cores",
    ):
        result.append(job.get())

    return result


def make_windows(pathname, args):
    splices = []
    idx = 0
    outfile = os.path.join(
        pathname, f"{os.path.splitext(os.path.basename(args.test))[0]}_split.fas"
    )
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
        query = os.path.join(query_dir, key + ".fasta")
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
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def parse_results(pathname):
    """Read % identity from Blast result files."""
    out, index = [], []

    for f in sorted(os.listdir(path=pathname)):
        index.append(f)
        curr = os.path.join(pathname, f)
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


def make_blast_db(pathname, filename):
    """Generate Blast database."""
    outfile = os.path.join(
        pathname, os.path.splitext(os.path.basename(filename))[0] + "-db"
    )
    cmd = (str(MAKEBLASTDB), "-dbtype", "nucl", "-in", filename, "-out", outfile)
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    return outfile


def main(args):
    logging.basicConfig(filename="fungani.log", level=logging.INFO)
    logger.info("creating temporary directories")
    if args.outdir is None:
        temp_dir = tempfile.TemporaryDirectory()
        fs_tmp = temp_dir.name
    else:
        fs_tmp = args.outdir
    fs_ani_queries = os.path.join(fs_tmp, "ani_q")
    os.makedirs(fs_ani_queries)
    fs_ani_blasts = os.path.join(fs_tmp, "ani_b")
    os.makedirs(fs_ani_blasts)

    logger.info("writing Blast db")
    reference_db = make_blast_db(fs_tmp, args.reference)

    logger.info("preparing sliced genome")
    fsplit = make_windows(fs_tmp, args)

    records = list(SeqIO.parse(fsplit, "fasta"))
    nsample = int(len(records) * args.percent / 100)
    records = random.sample(records, nsample)
    nc = int(len(records) / args.cpus)
    data = chunks(records, nc)
    p = multiprocessing.Pool(processes=args.cpus)

    logger.info("running Blast")
    result = run_async(
        blast,
        data,
        {"blast_dir": fs_ani_blasts, "query_dir": fs_ani_queries, "db": reference_db},
        args.cpus,
    )

    logger.info("parsing Blast")
    out, index = parse_results(fs_ani_blasts)
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

    if args.clean:
        logger.info("deleting temporary directories")
        shutil.rmtree(fs_tmp)
