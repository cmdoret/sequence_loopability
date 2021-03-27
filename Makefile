.PHONY: setup deps data

deps:
	python -m pip install -r requirements.txt

setup:
	python -m pip install -e .

data:
	dvc pull
