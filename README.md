# MSA Proteoform Profiler

**Goal**: Detect and annotate protein domains/isoforms from rna-seq multiple sequence alignments



### Dependencies

We use python 3.11 and [poetry](https://github.com/python-poetry/poetry) 1.6.1.
All package versions are managed by poetry.

### Running/Testing

Activate the virtual environment, then prepare by installing:
```
poetry shell
poetry install
```

**Running** the CLI:
```
poetry run msapp [options] <msa file>
```

**Running** the GUI:
```
poetry run msapp-gui
```

**Testing**:
```
poetry run pytest
```

Exit the venv:
```
poetry exit
```
