import cobra
import json


def json_load(path: str):
    """Loads the given JSON file and returns it as dictionary.

    Arguments
    ----------
    * path: str ~ The path of the JSON file
    """
    with open(path) as f:
        dictionary = json.load(f)
    return dictionary


dG0s = json_load("commodelpy/publication_runs/ecoli_models/publication_sbmls_and_dG0_jsons/ecolicore2double__dG0_ARB_no_h2o_no_pi__metabolites_periplasmic_cytosolic.json")
model = cobra.io.read_sbml_model("commodelpy/publication_runs/ecoli_models/publication_sbmls_and_dG0_jsons/ecolicore2double__dG0_ARB_no_h2o_no_pi__metabolites_periplasmic_cytosolic.xml")
model = cobra.io.read_sbml_model("commodelpy/publication_runs/ecoli_models/publication_sbmls_and_dG0_jsons/iML1515double__dG0_ARB_no_h2o_no_pi__metabolites_periplasmic_cytosolic.xml")
nans = 0
reactions = 0
for reaction in model.reactions:
    if reaction.id not in dG0s.keys():
        if reaction.id.endswith("ecoli1"):
            nans += 1
            reactions += 1
    else:
        if reaction.id.endswith("ecoli1"):
            reactions += 1

print("Number of NaNs:", nans)
print("Full number of internal reactions: ", reactions)
print("NaN fraction: ", (nans/reactions)*100, "%")
