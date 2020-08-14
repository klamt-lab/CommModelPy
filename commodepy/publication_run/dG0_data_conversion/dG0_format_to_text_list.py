from typing import Dict
from ...submodules.helper_general import json_write
import sys


# Separate from commodepy: Create dG0 dictionary for ecolicore2double and ecgsDouble
with open("./commodepy/publication_run/dG0_data_conversion/reaction_notes_iJO1366_with_extra_data.txt", "r") as f:
    dG0_lines = f.readlines()
dG0_lines = [x.replace("\n", "") for x in dG0_lines]

dG0_data_dict: Dict[str, str] = {}
for dG0_line in dG0_lines:
    reaction_name = "XXX" + dG0_line.split("    ")[0]
    reaction_name = reaction_name.replace("XXXR_", "")
    reaction_name = reaction_name.replace("XXX", "")

    if reaction_name in dG0_data_dict.keys():
        print("ERROR")
        sys.exit(-1)

    part_with_dG0 = dG0_line.split("deltaGR0;#;num;#;")[1]
    dG0 = part_with_dG0.split(";")[0]

    part_with_uncertainty = dG0_line.split("deltaGR0_Uncertainty;#;num;#;")[1]
    uncertainty_str = part_with_uncertainty.split(";")[0]
    uncertainty = float(uncertainty_str)

    dG0_data_dict[reaction_name] = dG0

output = "Reaction ID;dG0\n"
for key in dG0_data_dict.keys():
    output += key + ";" + dG0_data_dict[key] + "\n"

with open("./commodepy/publication_run/dG0_data_text_list/dG0.txt", "w") as f:
    f.write(output)


def run():
    """Function used for calling with "publication_run.py"""
    pass
