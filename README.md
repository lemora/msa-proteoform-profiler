# MSA Proteoform Profiler

**Goal**: Detect and annotate protein domains/isoforms from rna-seq multiple sequence alignments



### Dependencies

This software requires python >=3.9 and [poetry](https://github.com/python-poetry/poetry).
All package versions are managed by poetry.

### Running/Testing

Create the virtual environment with an appropriate python version (once):

```
poetry env use python3.12
```

Activate the virtual environment, then prepare by installing:
```
poetry shell
poetry install
```

**Running the CLI**:
```
msapp [options] <msa file>
```

Or, when not in venv: `poetry run msapp [options] <msa file>`.

**Running the GUI**:
```
msapp-gui
```

Or, when not in venv: `poetry run msapp-gui`.


**Testing**:
```
poetry run pytest
```

Exit the venv:
```
poetry exit
```
