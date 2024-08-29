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

- [Python ≥ 3.9](https://www.python.org/)
- [BLAST+ executables](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)
- [bedtools](https://bedtools.readthedocs.io/en/latest/index.html)

Python packages are managed using [Poetry](https://python-poetry.org/). External
packages will be installed automatically in a virtual environmenet when building
the project.

### Installing

FungANI should be available on [PyPI](https://pypi.org/) in the future. Assuming
you have a working installation of Python, use your package manager or simply
run the following instruction at a Python prompt:

    pip install --user FungANI

If you prefer to install the package from source, clone this repository

    git clone https://github.com/podo-gec/fungani.git

Or you can simply install it locally from GitHub as follows:

    pip install git+https://https://github.com/podo-gec/fungani

## Running the tests

There is a small test suite available in the tests directory. To run all the
tests, use pytest as follows:

    pytest tests/

Tests can be runned individually, see below.

### Sample tests

Some sample Fasta files are provided in order to check that everything works
fine.

    pytest tests/test_read_fasta.py
    pytest tests/test_read_collection.py

### Core tests

Blast and plotting capacilities are tested using a series of tests.

    pytest tests/test_blast_compat.py
    pytest tests/test_blast_results.py
    pytest tests/test_plot_pairwise.py
    pytest tests/test_plot_multiway.py

## Contributing

If you notice an unexpected result or if the application crashes, please fill an
issue and provide a bug report.

## Authors

- **Christophe Lalanne** ([chl@aliquote.org](mailto:chl@aliquote.org))
- **Philippe Silar** ([philippe.silar@u-paris.fr](mailto:philippe.silar@u-paris.fr))

## License

BSD-3
