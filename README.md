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

To use the R backend, please ensure that the `ggplot2` and `patchwork` packages
are installed system-wide. Alternatively, you could use `renv` or conda to
manage a shared virtual environment for both programs.

Python packages are managed using [Poetry]. External packages will be installed
automatically in a virtual environment when building the project. Currently,
the only external dependency is [fastaparser], which should be installed using
Poetry or pip install, e.g.

    python -m venv .venv # or source .venv/bin/activate if virtual env exists
    pip install -r requirements.txt

[Python ≥ 3.9]: https://www.python.org/
[BLAST+ executables]: https://blast.ncbi.nlm.nih.gov/doc/blast-help/
[R ≥ 4.0]: https://cran.r-project.org/
[conda]: https://docs.conda.io/en/latest/
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

Or you can simply install it locally from GitHub as follows:

    pip install git+https://https://github.com/podo-gec/fungani

[prerequisites]: #prerequisites

### Running the command-line application

Assuming you are at the root of the project and the virtual environment is
activated, simply run:

    python -m fungani.cli

Optionally, if R is installed on your OS, a graphical representation of the ANI
distribution will generated along raw results in your home user directory.

### Running the graphical application

Assuming you are at the root of the project and the virtual environment is
activated, simply run:

    python -m fungani.app

Optionally, if R is installed on your OS, a graphical representation of the ANI
distribution will generated along raw results in your home user directory.

The graphical application offers the exact same set of options as the
command-line application, see below.

![app](https://github.com/podo-gec/fungani/blob/master/assets/2024-08-30-14-39-38.png)

## Running the tests

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
