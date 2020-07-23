[![PyPI version](https://badge.fury.io/py/pyredcom-Paulocracy.svg)](https://badge.fury.io/py/pyredcom-Paulocracy)

# pyredcom

## Description

pyredcom is a collection of Python [dataclasses](https://docs.python.org/3/library/dataclasses.html) and associated functions which aim to help one in generating stoichiometric metabolic models of communities which consist of one or multiple species. All dataclasses and functions are based on [cobrapy](https://github.com/opencobra/cobrapy).

The underlying method for the definition of the generation of the community models and their MILP-based FBA analyses is the "RedCom" method as published in the following paper (Note: pyredcom's author is not directly associated with this publication, and pyredcom was not used for this paper):

* Koch, S., Kohrs, F., Lahmann, P., Bissinger, T., Wendschuh, S., Benndorf, D., ... & Klamt, S. (2019). RedCom: A strategy for reduced metabolic modeling of complex microbial communities and its application for analyzing experimental datasets from anaerobic digestion. *PLoS computational biology*, 15(2), e1006759. [doi:10.1371/journal.pcbi.1006759]( https://doi.org/10.1371/journal.pcbi.1006759)


## Documentation

A documentation of pyredcom's features and functions can be found in the "docs" subfolder of this repository. The documentation's starting point is "index.html". The whole documentation was generated using [pdoc3](https://github.com/pdoc3/).


## Installation procedure

### Option 1: As PyPI module

You can install pyredcom as Python module from [PyPI](https://pypi.org/project/autopacmen-Paulocracy/) using *pip*:

<pre>
pip install pyredcom-Paulocracy
</pre>

Afterwards, you can use pyredcom just as any other Python module using *import* in your Python session/script:
<pre>
import pyredcom
</pre>

You can also access pyredcom's documentation using Python's help function after importing pyredcom, e.g. for the whole module:
<pre>
help(pyredcom)
</pre>


### Option 2: Direct download

If you don't want to use pyredcom as PyPI module, you can also download this repository directly. pyredcom's main script file is the pyredcom Python script in the "pyredcom" subfolder.


## Exemplary usage

Exemplary usage of pyredcom can be found in the "puc" subfolder of the "pyredcom" folder. It contains the scripts for the generation
of the "shinjired" toy model, as well as community versions (i.e., a community of two times the respective model) of iJO1366 and
EColiCore2. See the "example"'s folder <span>README.md</span> for the references of these examples.


## Source code test range

A small test range for pyredcom can be found in the "test" subfolder of the "pyredcom" folder,
it is contained in the "test_pyredcom()" function of "test_pyredcom.py". The test system works e.g. with
[pytest](https://github.com/pytest-dev/pytest).


## License
pyredcom is free and open source, using the Apache License, Version 2
