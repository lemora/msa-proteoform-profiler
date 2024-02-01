# MSA Proteoform Profiler

The `MSAPP` is a tool for high-level analysis of RNA-Seq Multiple Sequence Alignments.
It assists in detecting and annotating proteoforms and domains.


### Dependencies

This software requires python >=3.9 and [poetry](https://github.com/python-poetry/poetry).
All package versions are managed by poetry.

### Running

Install poetry and python.
Create a virtual environment with an appropriate python version (once):

```
poetry env use python3.12
```

Activate the virtual environment, then prepare by installing:
```
poetry shell
poetry install
```

#### The Graphical User Interface

The GUI mode is currently the recommended way to use the MSAPP. It allows interactive analysis and fine-tuning of filtering parameters and similarity cutoff levels.

Running the GUI:
```
msapp-gui
```

Or, when not in venv: `poetry run msapp-gui`.


#### The Command-Line Interface

The CLI is currently a work in progress.

Running the CLI:
```
msapp [options] <msa file>
```

Or, when not in venv: `poetry run msapp [options] <msa file>`.

