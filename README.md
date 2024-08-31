# FungANI

A BLAST-based program for analyzing Average Nucleotide Identity (ANI) between
two fungal genomes, enables easy fungal species delimitation.

## Getting Started

These instructions will give you a copy of the project up and running on
your local machine for development and testing purposes. See deployment
for notes on deploying the project on a live system.

### Prerequisites

This software is built using Python and requires additional programs to be
installed separately.

- [Python ≥ 3.9]
- [BLAST+ executables]
- [R ≥ 4.0] (optional)

Standalone Blast programs can be installed on all platforms from the NCBI
website, or with [conda] on Linux and macOS systems.

To use the R backend, please ensure that the [ggplot2] and [patchwork] packages
are installed system-wide. Alternatively, you could use [renv] or conda to
manage a shared virtual environment for both programs.

Python packages are managed using [Poetry]. External packages will be installed
automatically in a virtual environment when building the project. Currently,
the only external dependency is [fastaparser], which should be installed using
Poetry or pip install.

[Python ≥ 3.9]: https://www.python.org/
[BLAST+ executables]: https://blast.ncbi.nlm.nih.gov/doc/blast-help/
[R ≥ 4.0]: https://cran.r-project.org/
[conda]: https://docs.conda.io/en/latest/
[ggplot2]: https://ggplot2.tidyverse.org/
[patchwaork]: https://patchwork.data-imaginist.com/index.html
[renv]: https://rstudio.github.io/renv/index.html
[Poetry]: https://python-poetry.org/
[fastaparser]: https://pypi.org/project/fastaparser/

### Installing

FungANI may be available on [PyPI](https://pypi.org/) in the future. In the
meantime, this package can be installed from from source by cloning the
repository from GitHub:

    git clone https://github.com/podo-gec/fungani.git

In this case, you will need to take care of creating and activating a virtual
environment, installing the relevant packages (see [prerequisites]), and
launching the application yourself (see below).

Example of use:

    git clone https://github.com/podo-gec/fungani.git
    python -m venv .venv
    source .venv/bin/activate
    pip install -r requirements.txt

Your application should be ready. You can deactivate the environment using the
corresponding command, anytime at your terminal prompt, when you are done. The
next time, you will only need to activate the virtual environment.

[prerequisites]: #prerequisites

### Running the command-line application

Assuming you are at the root of the project and the virtual environment is
activated, simply run:

    python -m fungani.cli REFERENCE TEST [-w 1000] [-g 500] [-j 20] [-c] [-o tmp]

In the above example, `REFERENCE` and `TEST` denote the path to the Fasta file
for the reference and test genomes. All other parameters are optional but you
can change window size (`-w` or `--window`), overlap (`-g` or `--overlap`),
number of cores (`-j` or `--cpus`), and output directory (`-o` or `--output`) to
store intermediate results. All intermediate results are cleaned up when the
application has finished its job.

If everything went fine, three files are written in your user home directory,
two CSV files that contain the % identity (on a 0-1 scale, i.e. 0.8 means 80%)
in the forward (`fungani_fwd.csv`) or reverse (`fungani_rev.csv`) mode. The
latter considers the `TEST` genome as the `REFERENCE` and all blasts are
performed against the `TEST` genome itself.

Optionally, if R is installed on your OS, a graphical representation of the ANI
distribution will be generated along raw results in your home user directory.

### Running the graphical application

Assuming you are at the root of the project and the virtual environment is
activated, simply run:

    python -m fungani.app

Optionally, if R is installed on your OS, a graphical representation of the ANI
distribution will be generated along raw results in your home user directory.

The graphical application offers the exact same set of options as the
command-line application, see above.

![app](https://github.com/podo-gec/fungani/blob/master/assets/2024-08-30-14-39-38.png)

## Performance

Expect application slowdowns on portable desktops. Performance issues are mainly
due to file read and write operations where Python does not really compete with
programs written in C/C++, Go or Rust. The two bottlenecks here are (1) splicing
the reference genome (since the [bedtools] suite is not available on Windows,
this is performed using builtin Python functions), and (2) writing intermediate
Fasta queries and Blast results as plain text files and reading them afterwards.
This, however, ensures a reasonably sized binary executable for Windows and
macOS users as we don't rely on [biopython] which would entail a large penalty
as it relies on [numpy]. Likewise, the graphical output is delegated to R as an
option, since including, e.g., [plotnine] would increase the binary size by
quite a large amount of Mb.

A more efficient command-line only version of this app, with memory-cached
operations for Linux and macOS users, will be available in a separate branch.

It should be noted, however, thta with the present version of the application,
it is possible to post-process the Blast results when the intermediate files are
kept (option `-c` or `--clean`). Two applications are possible with this setup:
(1) zero-hit regions can be analyzed individually, to highlight genomic
landscapr specific to each genomes; (2) using a window large enough to cover the
full sequence (e.g., with a databank of ITS) will result in a 'multi-blast'
setting whereby each sequence is blasted against a reference genome (in the
forward mode only).

Some benchmarks are shown below:

|  Processor                              | OS                   | No. cores (`-j`) | Time (HH:MM:SS)   |
| --------------------------------------- | -------------------- | ---------------- | ----------------- |
| Intel i7-10610U (8) @ 4.900GHz          | Ubuntu 24.04 LTS     | 4                | 00:32:32          |
| Intel Xeon Gold 6240R (96) @ 4.000GHz   | Ubuntu 22.04.4 LTS   | 20               | 00:12:36          |
| Intel Xeon Gold 6240R (96) @ 4.000GHz   | Ubuntu 22.04.4 LTS   | 40               | 00:08:02          |

In all cases CPU governor was set to "performance". Running in quick-mode (10%
of the genome) results in increased performance, especially on laptop (Intel
i7-10610U = 00:05:19 instead of 00:32:32).

Results from a [sample session] (whole genome and 10% sampling) are available.
Computations were performed on the genomes of _Neurospora crassa OR74A_ and
_Neurospora africana FGSC 1740_, available publicly on the NCBI databank.

[bedtools]: https://bedtools.readthedocs.io/en/latest/index.html
[biopython]: https://biopython.org/
[numpy]: https://numpy.org/
[plotnine]: https://plotnine.org/
[sample session]: https://github.com/podo-gec/fungani/blob/master/assets/sample_results

## Running the tests (dev-only)

There is a small test suite available in the tests directory. To run all the
tests, use pytest as follows:

    python -m pytest tests/

Alternatively, if you are using Poetry, since pytest is installed as a dev
dependency, you can simply run:

    poetry run pytest

Tests can be launched individually, see below.

### Sample tests

Some sample Fasta files are provided in order to check that everything works
fine.

    pytest tests/test_fasta_parser.py
    pytest tests/test_cli_parser.py

### Core tests

Blast and plotting capabilities are also tested using a series of tests.

## Contributing

If you notice an unexpected result or if the application crashes, please fill an
issue and provide a bug report, or contact the authors by email (see below).

## Authors

- **Christophe Lalanne** ([chl@aliquote.org])
- **Philippe Silar** ([philippe.silar@u-paris.fr])

[chl@aliquote.org]: mailto:chl@aliquote.org
[philippe.silar@u-paris.fr]: mailto:philippe.silar@u-paris.fr

## License

BSD-3
