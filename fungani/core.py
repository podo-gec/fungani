import csv
import datetime
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
import time

import fastaparser

logger = logging.getLogger(__name__)
BLASTN = shutil.which("blastn")
MAKEBLASTDB = shutil.which("makeblastdb")


def deserialize_fasta(sequences):
    records = []
    for seq in sequences:
        records.append([seq.sequence, seq.header.strip(">")])
    return records


def run_async(func, arglist, kwds, cpus):
    pool = multiprocessing.Pool(processes=cpus)
    jobs = [
        pool.apply_async(
            func=func,
            args=(arg,),
            kwds=kwds,
        )
        for arg in arglist
    ]

    pool.close()
    # pool.join()
    results = []
    for job in jobs:
        results.append(job.get())

    return results


def make_windows(pathname, args):
    splices = []
    idx = 0
    outfile = os.path.join(
        pathname, f"{os.path.splitext(os.path.basename(args.test))[0]}_split.fas"
    )

    with open(args.test, "r") as file:
        records = list(fastaparser.Reader(file, parse_method="quick"))

    for record in records:
        idx += 1
        end = len(record.sequence) - args.size + 1
        for start in range(0, end, args.overlap):
            id = ":".join((str(idx), "-".join((str(start), str(start + args.size)))))
            seq = record.sequence[start : start + args.size]
            splices.append(fastaparser.FastaSequence(seq, id))

    with open(outfile, "w") as file:
        writer = fastaparser.Writer(file)
        for splice in splices:
            writer.writefasta(splice)

    return outfile, len(splices)


def chunks(lst, n):
    if n == 0:
        n = 1
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


def blast(record, blast_dir, query_dir, db):
    for r in record:
        key = str(r[1])
        query = os.path.join(query_dir, key + ".fasta")
        with open(query, "w") as file:
            writer = fastaparser.Writer(file)
            writer.writefasta(fastaparser.FastaSequence(r[0], r[1]))
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


def main(args, start_time=None):
    logging.basicConfig(
        filename=os.path.join(os.path.expanduser("~"), "fungani.log"),
        level=logging.INFO,
        filemode="w+" if args.mode == "fwd" else "w",
    )

    if args.mode == "fwd":
        logger.info("========== Starting new process ==========")
    logger.info(f"==========       mode: {args.mode}      ==========")

    logger.info("Creating temporary directories")
    if args.outdir is None:
        temp_dir = tempfile.TemporaryDirectory()
        fs_tmp = temp_dir.name
    else:
        if os.path.isdir(args.outdir):
            if os.listdir(args.outdir):
                sys.exit(f"Directory {args.outdir} is not empty")
        fs_tmp = args.outdir
    fs_ani_queries = os.path.join(fs_tmp, "ani_q")
    os.makedirs(fs_ani_queries)
    fs_ani_blasts = os.path.join(fs_tmp, "ani_b")
    os.makedirs(fs_ani_blasts)

    if args.size <= 15:
        sys.exit(f"Window size {args.size} too small")

    if args.cpus > multiprocessing.cpu_count() - 1:
        args.cpus = multiprocessing.cpu_count() - 1
        logger.warn(f"No. cores too large, now {args.cpus}")

    if 0 < args.percent <= 1:
        args.percent *= 100
        logger.warn(f"Genome sampling < 1, now {args.percent}")

    logger.info("Writing Blast database")
    reference_db = make_blast_db(fs_tmp, args.reference)

    logger.info("Preparing sliced genome")
    fsplit, count = make_windows(fs_tmp, args)
    logger.info(f">>> {count} sequences in spliced genome")

    with open(fsplit, "r") as file:
        records = list(fastaparser.Reader(file, parse_method="quick"))
    records = deserialize_fasta(records)
    nsample = int(len(records) * args.percent / 100)
    records = random.sample(records, nsample)
    nc = int(len(records) / args.cpus)
    data = chunks(records, nc)

    logger.info("Running Blast")
    result = run_async(
        blast,
        data,
        {"blast_dir": fs_ani_blasts, "query_dir": fs_ani_queries, "db": reference_db},
        args.cpus,
    )
    logger.info(f">>> {len(result)} jobs terminated successfully")

    logger.info("Parsing Blast")
    out, index = parse_results(fs_ani_blasts)
    logger.info(f">>> {len(out)} Blast results analysed")
    outfile = os.path.join(os.path.expanduser("~"), f"fungani_{args.mode}.csv")
    rows = zip(index, out)
    with open(outfile, "w") as f:
        writer = csv.writer(f)
        for row in rows:
            writer.writerow(row)
    logger.info(">>> results saved in user home directory")

    if args.clean:
        logger.info("Deleting temporary directories")
        shutil.rmtree(fs_tmp)
        logger.info(f">>> '{fs_tmp}' deleted successfully")

    if args.mode == "rev":
        logger.info("=========== Process completed ============")
        if start_time is not None:
            toc = time.time()
            elapsed_time = datetime.datetime.strftime(
                datetime.datetime.fromtimestamp(toc - start_time, datetime.UTC),
                "%H:%M:%S",
            )
            logger.info(f"Elapsed time: {elapsed_time}")
            # Simulate a success return value for the Tk app
            return True
