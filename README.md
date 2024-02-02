# MSA Proteoform Profiler

The `MSAPP` is a tool for high-level analysis of RNA-Seq multiple sequence alignments.
It assists in detecting and annotating proteoforms and domains.

It can either be run as with graphical user interface (recommended), or as a command-line application.


### Dependencies

This software requires python >= 3.9 and [poetry](https://github.com/python-poetry/poetry).
All package versions are managed by poetry.

### Running

Install poetry. Create a virtual environment with an appropriate python version (once):

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

The CLI offers similar tuning/visualization parameters as the GUI, but it is less interactive.

Running the CLI:
```
msapp [options] <msa file>
```

Or, when not in venv: `poetry run msapp [options] <msa file>`.

