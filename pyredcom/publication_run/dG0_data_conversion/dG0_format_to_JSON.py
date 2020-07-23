from typing import Dict
from ...submodules.helper_general import json_write
import sys


# Separate from pyredcom: Create dG0 dictionary for ecolicore2double and ecgsDouble
with open("./pyredcom/publication_run/dG0_data_conversion/reaction_notes_iJO1366_with_extra_data.txt", "r") as f:
    dG0_lines = f.readlines()
dG0_lines = [x.replace("\n", "") for x in dG0_lines]

dG0_data_dict: Dict[str, Dict[str, float]] = {}
for dG0_line in dG0_lines:
    reaction_name = "XXX" + dG0_line.split("    ")[0]
    reaction_name = reaction_name.replace("XXXR_", "")
    reaction_name = reaction_name.replace("XXX", "")

    if reaction_name in dG0_data_dict.keys():
        print("ERROR")
        sys.exit(-1)

    part_with_dG0 = dG0_line.split("deltaGR0;#;num;#;")[1]
    dG0_str = part_with_dG0.split(";")[0]
    if dG0_str == "NaN":
        continue
    dG0_data_dict[reaction_name+"_ecoli1"] = {}
    dG0_data_dict[reaction_name+"_ecoli2"] = {}
    dG0 = float(dG0_str)

    part_with_uncertainty = dG0_line.split("deltaGR0_Uncertainty;#;num;#;")[1]
    uncertainty_str = part_with_uncertainty.split(";")[0]
    uncertainty = float(uncertainty_str)

    dG0_data_dict[reaction_name+"_ecoli1"]["dG0"] = dG0
    dG0_data_dict[reaction_name+"_ecoli2"]["dG0"] = dG0
    dG0_data_dict[reaction_name+"_ecoli1"]["uncertainty"] = uncertainty
    dG0_data_dict[reaction_name+"_ecoli2"]["uncertainty"] = uncertainty

json_write("./pyredcom/publication_run/dG0_data_JSON/dG0.json", dG0_data_dict)


def run():
    """Function used for calling with "publication_run.py"""
    pass
