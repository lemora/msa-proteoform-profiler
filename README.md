The MSA Proteoform Profiler
---

![Python3](https://img.shields.io/badge/Language-Python3-steelblue)
![OS](https://img.shields.io/badge/OS-_Windows_|_Mac_|_Linux-steelblue)
![License](https://img.shields.io/badge/License-GPLv3-steelblue)

Overview
---

The `MSAPP` is a tool for high-level analysis of RNA-Seq multiple sequence alignments.
It assists in detecting and annotating proteoforms and domains.

It can either be run with a graphical user interface (recommended), or as a command-line tool.


Dependencies
---

This software requires python >= 3.9 and [poetry](https://github.com/python-poetry/poetry).
All package versions are managed by poetry.

Running
---

Install poetry. Create a virtual environment with an appropriate python version (once):
`poetry env use python3.12`.

Activate the virtual environment: `poetry shell`. Then prepare by installing: `poetry install`.

### The Graphical User Interface

The GUI mode is currently the recommended way to use the MSAPP. It allows interactive analysis and fine-tuning of filtering parameters and similarity cutoff levels.

Running the GUI:
```
msapp-gui
```

Or, when not in venv: `poetry run msapp-gui`.


### The Command-Line Interface

The CLI offers similar tuning/visualization parameters as the GUI, but it is less interactive.

Running the CLI:
```
msapp [options] <msa file>
```

Or, when not in venv: `poetry run msapp [options] <msa file>`.


License
---

This software is licensed under GPLv3.0 or later.