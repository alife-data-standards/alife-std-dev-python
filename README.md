# ALife Data Standards - Python Development Utilities

[![CI](https://github.com/alife-data-standards/alife-std-dev-python/actions/workflows/CI.yaml/badge.svg)](https://github.com/alife-data-standards/alife-std-dev-python/actions/workflows/CI.yaml)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/0305e8a814e04d4395c25a70b2908651)](https://www.codacy.com/gh/alife-data-standards/alife-std-dev-python/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=alife-data-standards/alife-std-dev-python&amp;utm_campaign=Badge_Grade)
[![codecov](https://codecov.io/gh/alife-data-standards/alife-std-dev-python/branch/master/graph/badge.svg?token=FGMQICJ2SK)](https://codecov.io/gh/alife-data-standards/alife-std-dev-python)
![PyPI](https://img.shields.io/pypi/v/ALifeStdDev?color=blue)

This is the repository for the ALifeStdDev Python package, which contains Python
development utilities for working with [standardized](https://github.com/alife-data-standards/alife-data-standards)
ALife data.

# Installation Instructions

ALifeStdDev can be installed using pip:

```
pip install ALifeStdDev
```

# Usage Instructions

To load a single submodule,

```python3
from ALifeStdDev import phylogeny as asd_phylo
asd_phylo.load_phylogeny_to_pandas_df("myfile.csv")
```

To load the library as a flat namespace,

```python3
from ALifeStdDev import ALifeStdDev as asd
asd.load_phylogeny_to_pandas_df("myfile.csv")
```
