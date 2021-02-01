#!/usr/bin/env python3
#
# Copyright 2020-2021 PSB
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Generation script for the community model 'ecolicore2double'.

ecolicore2double models a community of two EcoliCore2 'organisms'. Its intention and usage is explained in
more detail in commmodelpy's publication.
"""
# IMPORT SECTION
# External modules
import cobra
import copy
from typing import Dict
# Internal modules
from commmodelpy.commmodelpy import Community, SingleModel, generate_community_model_with_no_growth
from commmodelpy.submodules.helper_general import json_write, json_load


# ecolicore double model
# with internal exchanges for everything in the periplasm
ecoli_model = cobra.io.read_sbml_model(
    "./publication_runs/ecoli_models/original_sbml_models_in_cleaned_form/ecolicore2_loaded_and_saved_by_cobrapy_cleaned.xml")

# Define essential and periplasmic metabolites (essential metabolites
# were found by FBAs beforehand)
essential_in_metabolites = [
    "co2_p",
    "fe2_p",
    "nh4_p",
    "cu2_p",
    "glc__D_p",
    "o2_p",
    "pi_p",
    "so4_p",
    "h_p",
    "h2o_p",
    # New ones
    "fe3_p",
    "mn2_p",
    "zn2_p",
    "mg2_p",
    "ca2_p",
    "ni2_p",
    "cobalt2_p",
    "mobd_p",
    "k_p",
    "cl_p",
]
essential_out_metabolites = [
    "co2_p",
    "ac_p",
    "h_p",
    "h2o_p",
    "etoh_p",
    "succ_p",
    "lac__D_p",
    "for_p",
    # "Demand" metabolites of DM_R reaction
    "4crsol_c",
    "5drib_c",
    "amob_c",
    "mththf_c",
    # Quasi "demand metabolite"
    "meoh_p",
]
all_essential_metabolites = list(
    set(essential_in_metabolites+essential_out_metabolites))
essential_metabolite_mapping: Dict[str, str] = {}
for met in essential_in_metabolites+essential_out_metabolites:
    essential_metabolite_mapping[met] = (
        met+"\b").replace("_p\b", "").replace("_c\b", "")

periplasmic_metabolites = [
    x for x in ecoli_model.metabolites if x.id.endswith("_p")]
periplasmic_metabolites_cut = [(x.id+"\b").replace("_p\b", "")
                               for x in periplasmic_metabolites]
periplasmic_metabolites += [x for x in ecoli_model.metabolites
                            if x.id.endswith("_c") and (((x.id+"\b").replace("_c\b", "") not in periplasmic_metabolites_cut))]
periplasmic_metabolites = [x for x in periplasmic_metabolites
                           if (((x.id+"\b").replace("_c\b", "_p") not in all_essential_metabolites)) and
                              (((x.id+"\b").replace("_p\b", "_p") not in all_essential_metabolites))]
all_input_metabolite_ids = []
all_output_metabolite_ids = []
all_inout_metabolite_ids_mapping: Dict[str, str] = {}
for periplasmic_metabolite in periplasmic_metabolites:
    metabolite_id = periplasmic_metabolite.id
    all_input_metabolite_ids.append(metabolite_id)
    all_output_metabolite_ids.append(metabolite_id)
    all_inout_metabolite_ids_mapping[metabolite_id] = (
        metabolite_id+"\b").replace("_p\b", "").replace("_c\b", "")

combined_input_metabolite_ids = list(set(
    all_input_metabolite_ids + essential_in_metabolites + essential_out_metabolites))
combined_output_metabolite_ids = list(set(
    all_output_metabolite_ids + essential_out_metabolites + essential_in_metabolites))
combined_metabolite_ids_mapping = {
    **all_inout_metabolite_ids_mapping, **essential_metabolite_mapping}

# Define commmodelpy SingleModel
ecoli_1 = SingleModel(
    cobra_model=ecoli_model,
    species_abbreviation="ecoli1",
    objective_reaction_id="Ec_biomass_iJO1366_core_53p95M",
    exchange_reaction_id_prefix="EX_",
    input_metabolite_ids=combined_input_metabolite_ids,
    output_metabolite_ids=combined_output_metabolite_ids,
    model_metabolite_to_exchange_id_mapping=combined_metabolite_ids_mapping,
)
# Double it :D
ecoli_2 = copy.deepcopy(ecoli_1)
ecoli_2.species_abbreviation = "ecoli2"

# Create community model :D
print("===\nGeneration of community model...")
potential_product_metabolites = [x for x in periplasmic_metabolites]
potential_product_metabolite_ids = [
    (x.id+"\b").replace("_p\b", "").replace("_c\b", "") for x in potential_product_metabolites]

community = Community(
    single_models=[ecoli_1, ecoli_2],
    exchange_compartment_id="exchg",
    exchange_reaction_id_prefix="EX_C_",
    input_metabolite_ids=[(x+"\b").replace("_p\b", "")
                          for x in essential_in_metabolites],
    output_metabolite_ids=[(x+"\b").replace("_p\b", "").replace("_c\b", "").replace("\b", "")
                           for x in essential_out_metabolites+potential_product_metabolite_ids]
)
community_model = generate_community_model_with_no_growth(
    community, {"ecoli1": 0.5, "ecoli2": 0.5})

# Separate from commmodelpy: Create dG0 dictionary for ecolicore2double
dG0_data_dict = json_load(
    "./publication_runs/ecoli_models/model_reaction_id_to_dG0_mapping_jsons/reaction_id_to_dG0_mapping_for_ecolicore.json")


thermodynamically_excluded_metabolites = ["h2o", "h"]
for reaction in community_model.reactions:
    if reaction.id.startswith("EXCHG_"):
        exclude_reaction = False
        for excluded_metabolite in thermodynamically_excluded_metabolites:
            if reaction.id.endswith("_to_"+excluded_metabolite):
                exclude_reaction = True
        if not exclude_reaction:
            dG0_data_dict[reaction.id] = {}
            dG0_data_dict[reaction.id]["dG0"] = 0
            dG0_data_dict[reaction.id]["uncertainty"] = 0
        is_essential = False
        for essential_metabolite in all_essential_metabolites:
            if reaction.id.endswith("_to_"+(essential_metabolite+"\b").replace("_p\b", "").replace("_c\b", "")):
                is_essential = True
                break
        if is_essential:
            continue
        reaction.lower_bound = 0
        reaction.upper_bound = 0
    elif reaction.id.startswith("EX_C"):
        is_essential = False
        for essential_metabolite in all_essential_metabolites:
            if reaction.id == ("EX_C_"+(essential_metabolite+"\b").replace("_p\b", "_exchg").replace("_c\b", "_exchg")):
                is_essential = True
                break
        if is_essential:
            print(reaction.id, "is recognized as essential")
            continue
        reaction.lower_bound = 0
        reaction.upper_bound = 0

print("Done!")

# Store model as SBML :D
cobra.io.write_sbml_model(
    community_model, "./publication_runs/ecoli_models/publication_sbmls_and_dG0_jsons/ecolicore2double_model.xml")
json_write("./publication_runs/ecoli_models/publication_sbmls_and_dG0_jsons/ecolicore2double_dG0.json", dG0_data_dict)
