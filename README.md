# sequence_loopability
Investigating the sequence features underlying chromatin loops.

## Project structure:

* src: Stores the code for the main pipeline, data extraction, processing, training, evaluation...
* notebooks: Exploratory analyses in the form of jupyter notebooks.
* seqloops: Boilerplate code and utilities meant to be imported as a python package.
* scripts: various scripts that were used to generate the input data.


## Setup:

To make `seqloops` importable in python scripts and notebooks, you can run: `make setup`.

All input and output data are managed via dvc. They can be imported as follows:
```bash
pip install dvc[gdrive]
dvc pull
```

## Workflow

Code changes are managed via `git`, data changes are managed via `dvc`.
