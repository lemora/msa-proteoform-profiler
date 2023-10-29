# MSA Proteoform Profiler

**Goal**: Detect and annotate protein domains/isoforms from rna-seq multiple sequence alignments



### Dependencies

We use python 3.11 and [poetry](https://github.com/python-poetry/poetry) 1.6.1.
All package versions are managed by poetry. It is advised to create a virtual environment:
```
poetry env use python3
poetry env info
```

### Testing

Activate the viurtual environment with `poetry shell`, then:
```
poetry install
poetry run pytest
```
Exit the venv with `poetry exit`.
