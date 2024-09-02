.PHONY: tests
.DEFAULT_GOAL := tests

check:
	ruff check fungani

format:
	black fungani

tests:
	python -m pytest -vv

build:
	poetry build
