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
are installed system-wide. Alternatively, you could use [renv] or [conda] to
manage a shared virtual environment for both programs.

Python packages are managed using [Poetry]. External packages will be installed
automatically in a virtual environment when building the project. Currently,
the only external dependency is [fqfa], which should be installed using
Poetry or pip install.

[Python ≥ 3.9]: https://www.python.org/
[BLAST+ executables]: https://blast.ncbi.nlm.nih.gov/doc/blast-help/
[R ≥ 4.0]: https://cran.r-project.org/
[conda]: https://docs.conda.io/en/latest/
[ggplot2]: https://ggplot2.tidyverse.org/
[patchwaork]: https://patchwork.data-imaginist.com/index.html
[renv]: https://rstudio.github.io/renv/index.html
[conda]: https://www.anaconda.com/download/
[Poetry]: https://python-poetry.org/
[fqfa]: https://pypi.org/project/fqfa/

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
    cd fungani
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

    python -m fungani.cli -h

This should print the help message and show options with default values:

    usage: fungani [-h] [-t THRESHOLD] [-p PERCENT] [-w SIZE] [-g OVERLAP] [-j CPUS] [-o OUTDIR] [-u] [-c] reference test

    Compute ANI using Blast

    positional arguments:
      reference                            Reference genome
      test                                 Test genome

    options:
      -h, --help                           show this help message and exit
      -t THRESHOLD, --threshold THRESHOLD  ANI threshold (default: 80)
      -p PERCENT, --percent PERCENT        Genome fraction (default: 10)
      -w SIZE, --size SIZE                 Window size (default: 1000)
      -g OVERLAP, --overlap OVERLAP        Window overlap (default: 500)
      -j CPUS, --cpus CPUS                 Number of CPU cores (default: 4)
      -o OUTDIR, --output OUTDIR           Output directory (default: None)
      -u, --onepass                        Only in one direction only (default: False)
      -c, --clean                          Clean intermediate files (default: False)

To run the program on a REFERENCE and TEST genome, use the following command:

    python -m fungani.cli REFERENCE TEST [-t 80] [-p 10] [-w 1000] [-g 500] [-j 20] [-c] [-o tmp]

In the above example, `REFERENCE` and `TEST` denote the path to the Fasta file
for the reference and test genomes. All other parameters are optional and can
safely be omitted but you can change ANI threshold (`-t` or `--threshold`),
genome fraction (`-p` or `--percent`), window size (`-w` or `--window`), overlap
(`-g` or `--overlap`), number of cores (`-j` or `--cpus`), direction (`-u` or
`--onepass`), cleaning of intermediate results (̀`c` or `--clean`), and output
directory (`-o` or `--output`) to store intermediate results. All intermediate
results are cleaned up when the application has finished its job.

If everything went fine, three files are written in your user home directory,
two CSV files that contain the % identity (on a 0-1 scale, i.e. 0.8 means 80%)
in the forward (`fungani_fwd.csv`) and reverse (`fungani_rev.csv`) direction if
you didn't activate the `-u` or `--onepass` option. The latter considers the
`TEST` genome as the `REFERENCE` and all blasts are performed against the `TEST`
genome itself.

Optionally, if R is installed on your OS, a graphical representation of the ANI
distribution will be generated along raw results in your home user directory.

### Running the graphical application

Assuming you are at the root of the project and the virtual environment is
activated, simply run:

    python -m fungani

Optionally, if R is installed on your OS, a graphical representation of the ANI
distribution will be generated along raw results in your home user directory.

The graphical application offers the exact same set of options as the
command-line application, see above.

![app](https://github.com/podo-gec/fungani/blob/master/assets/2024-09-03-15-07-27.png)

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

It should be noted, however, that with the present version of the application,
it is possible to post-process the Blast results when the intermediate files are
kept (option `-c` or `--clean`). Two applications are possible with this setup:
(1) zero-hit regions can be analyzed individually, to highlight genomic
landscape specific to each genomes; (2) using a window large enough to cover the
full sequence (e.g., with a databank of ITS) will result in a 'multi-blast'
setting whereby each sequence is blasted against a reference genome (in the
forward mode only).

Some benchmarks are shown below. In all cases CPU governor was set to
"performance". Computations were performed on the genomes of _Neurospora crassa_
OR74A and _Neurospora africana_ FGSC 1740, available publicly on the NCBI
databank.

Whole-genome analysis:

|  Processor                            | OS                 | No. cores | Time (HH:MM:SS) |
| ------------------------------------- | ------------------ | --------- | --------------- |
| Intel i7-10610U (8) @ 4.900GHz        | Ubuntu 24.04 LTS   | 4         | 00:22:32        |
| Intel Xeon E5-2630 v3 (32) @ 3.200GHz | Ubuntu 20.04.6 LTS | 20        | 00:07:53        |
| Intel Xeon E5-2630 v3 (32) @ 3.200GHz | Ubuntu 20.04.6 LTS | 30        | 00:07:05        |
| Intel Xeon Gold 6240R (96) @ 4.000GHz | Ubuntu 22.04.4 LTS | 20        | 00:07:56        |
| Intel Xeon Gold 6240R (96) @ 4.000GHz | Ubuntu 22.04.4 LTS | 40        | 00:04:01        |

Quick mode analysis (10% genome, `-p 10`):

|  Processor                            | OS                 | No. cores | Time (HH:MM:SS) |
| ------------------------------------- | ------------------ | --------- | --------------- |
| Intel i7-10610U (8) @ 4.900GHz        | Ubuntu 24.04 LTS   | 4         | 00:02:44        |
| Intel Xeon E5-2630 v3 (32) @ 3.200GHz | Ubuntu 20.04.6 LTS | 20        | 00:01:05        |
| Intel Xeon E5-2630 v3 (32) @ 3.200GHz | Ubuntu 20.04.6 LTS | 30        | 00:00:58        |
| Intel Xeon Gold 6240R (96) @ 4.000GHz | Ubuntu 22.04.4 LTS | 20        | 00:01:05        |
| Intel Xeon Gold 6240R (96) @ 4.000GHz | Ubuntu 22.04.4 LTS | 40        | 00:00:41        |

Results from a [sample session] (whole genome and 10% sampling) are available.

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

Note that tests can be launched individually.

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
