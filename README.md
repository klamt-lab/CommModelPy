[![PyPI version](https://badge.fury.io/py/commmodelpy.svg)](https://badge.fury.io/py/commmodelpy)

# CommModelPy

## Description

CommModelPy is a collection of Python [dataclasses](https://docs.python.org/3/library/dataclasses.html) and associated functions which aim to help one in generating stoichiometric metabolic models of communities which consist of one or multiple species. All dataclasses and functions are based on [cobrapy](https://github.com/opencobra/cobrapy).

The underlying methods for the generation and analysis of the community models is explained in more detail in CommModelPy's source code documentation.

An exemplary usage of CommModelPy is given in its publiation (Bekiaris & Klamt, in submission).

## Installation procedure

### Option 1: As PyPI module

You can install CommModelPy as Python module from [PyPI](https://pypi.org/project/commmodelpy/) using *pip*:

<pre>
pip install commmodelpy
</pre>

Afterwards, you can use CommModelPy just as any other Python module using *import* in your Python session/script:
<pre>
import commmodelpy
</pre>
In order to import CommModelPy's main script with all major dataclasses and functions and which is called "commmodelpy.py", you can import it using
<pre>
import commmodelpy.commmodelpy
</pre>

### Option 2: Direct download

If you don't want to use CommModelPy as PyPI module, you can also download this repository directly. The main script file is the commmodelpy.py Python script in the "commmodelpy" subfolder.

## Repository structure

* The actual commmodelpy pip package can be found in the "commmodelpy" subfolder, where "commmodelpy.py" contains all relevant functions and dataclasses.
* All Python scripts which were used in CommModelPy's publication, which use CommModelPy's function for community models without defined growth and with fixed species ratios, can be found in the "publication_runs" subfolder, which is in the "commmodelpy" subfolder. The scripts in the local subfolder "toy_model" contain the script for the generation of the toy model shown in the publication. The scripts in the local subfolder "ecoli_models" contain the generation of dG0 data using the [Equilibrator API](https://gitlab.com/equilibrator/equilibrator-api) as well as the CommModelPy-assisted generation of iML1515 and EcoliCore2 single-species community models. A complete call of all E. coli model scripts in the right order is given by the "execute_publication_ecoli_model_scripts.py" script in the main folder, a call of the toy model scripts is given by the "execute_publication_toy_model_script.py" script in the main folder.
* An exemplary usage of CommModelPy with its function with a defined fixed growth rate and free species ratios can be found in the "balanced_growth_example" subfolder. A call of the relevant script is given in the "execute_balanced_growth_example.py" script in the main folder.

## Documentation

A documentation of CommModelPy's features and functions can be found in the "docs/commmodelpy/" subfolder of this repository. The documentation's starting point is "index.html". The whole documentation was generated using [pdoc3](https://github.com/pdoc3/).

You can also access CommModelPy's documentation using Python's help function after importing CommModelPy, e.g. for the whole module:
<pre>
help(commmodelpy)
</pre>

## Publication

CommModelPy is published in the following publication:

* Bekiaris & Klamt, 2021, *in submission*

## License

CommModelPy is free and open source, using the Apache License, Version 2
