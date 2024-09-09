.PHONY: tests
.DEFAULT_GOAL := run

archive:
	tar cf fungani.tgz fungani/* poetry.lock requirements.txt pyproject.toml

check:
	ruff check fungani

format:
	black fungani

test:
	python -m pytest -vv

run:
	python -m fungani

build:
	poetry build
