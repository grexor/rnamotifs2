#!/bin/bash
# https://packaging.python.org/en/latest/tutorials/packaging-projects/
rm -f dist/*
python -m build
pip install .
