# sequence_loopability
Investigating the sequence features underlying chromatin loops.

## Project structure:

* src: Stores the code for the main pipeline, data extraction, processing, training, evaluation...
* notebooks: Exploratory analyses in the form of jupyter notebooks.
* seqloops: Boilerplate code and utilities meant to be imported as a python package.
* scripts: various scripts that were used to generate the input data.


## Setup:

All dependencies can be installed using:
```bash
make deps
```

To make `seqloops` importable in python scripts and notebooks, you can run: `make setup`.

All input and output data are managed via dvc. They can be imported as follows:
```bash
pip install dvc[gdrive]
dvc pull
```

## Workflow

Code changes are managed via `git`. Data changes are managed via `dvc`, which is connected to a google drive folder.
When modifying or adding new datafiles in the `data` folder, the modifications must be uploaded to the dvc server.
The updated small tracker file (`.dvc`) must be commited to git to keep track of changes.
The standard process is as follows:
```bash
dvc add data
dvc push
git add data.dvc
git commit -m 'added new file'
git push
```
