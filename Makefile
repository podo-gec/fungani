.PHONY: tests
.DEFAULT_GOAL := run

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
