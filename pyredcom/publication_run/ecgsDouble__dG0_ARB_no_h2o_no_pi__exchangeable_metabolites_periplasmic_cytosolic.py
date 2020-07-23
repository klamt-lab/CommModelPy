#!/usr/bin/env python3
#
# Copyright 2020 PSB
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
"""Generation script for the community model 'ecolicore2double'."""

# IMPORT SECTION
# External modules
import cobra
import copy
from typing import Dict
# Internal modules
from ..pyredcom import Community, SingleModel, generate_community_cobra_model, redcom_fba
from ..submodules.helper_general import json_write, json_load


# ecolicore double model
# with internal exchanges for everything in the periplasm
ecoli_model = cobra.io.read_sbml_model("./pyredcom/publication_run/original_sbml_models_loaded_and_saved_by_cobrapy/ecgs_loaded_and_saved_by_cobrapy.xml")
# ecoli_model.solver = "cplex"

# Solve wrong reaction name parts
for reaction in ecoli_model.reactions:
    if "_DASH_" in reaction.id:
        reaction.id = reaction.id.replace("_DASH_", "__")
    if "_LPAREN_" in reaction.id:
        reaction.id = reaction.id.replace("_LPAREN_", "_")
    if "_RPAREN_" in reaction.id:
        reaction.id = reaction.id.replace("_RPAREN_", "")

# Delete boundary-condition-free exchange metabolites
metabolite_ids = [x.id for x in ecoli_model.metabolites]
for metabolite_id in metabolite_ids:
    if metabolite_id.endswith("_ex"):
        metabolite_instance = ecoli_model.metabolites.get_by_id(metabolite_id)
        ecoli_model.remove_metabolites([metabolite_instance])

reaction_ids = [x.id for x in ecoli_model.reactions]
for reaction_id in reaction_ids:
    reaction = ecoli_model.reactions.get_by_id(reaction_id)
    if reaction.metabolites == {}:
        ecoli_model.remove_reactions([reaction])

# Rename wrong metabolite name parts
for metabolite in ecoli_model.metabolites:
    if "_DASH_" in metabolite.id:
        metabolite.id = metabolite.id.replace("_DASH_", "__")
    if "_LPAREN_" in metabolite.id:
        metabolite.id = metabolite.id.replace("_LPAREN_", "_")
    if "_RPAREN_" in metabolite.id:
        metabolite.id = metabolite.id.replace("_RPAREN_", "")

# Deactivate all C sources except of glucose
ecoli_model.reactions.EX_succ_e.lower_bound = 0
ecoli_model.reactions.EX_glyc_e.lower_bound = 0
ecoli_model.reactions.EX_ac_e.lower_bound = 0

# Test renamed model with FBA
with ecoli_model:
    print("Single model FBA solution:")
    fba_solution = ecoli_model.optimize()
    print(ecoli_model.summary())
    for reaction in ecoli_model.reactions:
        if not reaction.id.startswith("EX_"):
            continue
        if fba_solution.fluxes[reaction.id] != 0:
            print(f"{reaction.id}: {fba_solution.fluxes[reaction.id] }")
    print("~~~")


# Define essential and periplasmic metabolites (essential metabolites
# were found by FBAs beforehand)
essential_in_metabolites = [
    "fe2_e",
    "nh4_e",
    "cu2_e",
    "glc__D_e",
    "o2_e",
    "pi_e",
    "so4_e",
    "h_e",
    "h2o_e",
]
essential_out_metabolites = [
    "co2_e",
    "meoh_e",
    "ac_e",
    "h_e",
    "h2o_e",
    "etoh_e",
    "succ_e",
    "lac__D_e",
    "for_e",
]
all_essential_metabolites = list(set(essential_in_metabolites+essential_out_metabolites))
essential_metabolite_mapping: Dict[str, str] = {}
for met in essential_in_metabolites+essential_out_metabolites:
    essential_metabolite_mapping[met] = (met+"\b").replace("_e\b", "")

periplasmic_metabolites = [x for x in ecoli_model.metabolites if x.id.endswith("_p")]
periplasmic_metabolites_cut = [(x.id+"\b").replace("_p\b", "") for x in periplasmic_metabolites]
periplasmic_metabolites += [x for x in ecoli_model.metabolites
                            if x.id.endswith("_c") and (((x.id+"\b").replace("_c\b", "") not in periplasmic_metabolites_cut))]
periplasmic_metabolites = [x for x in periplasmic_metabolites
                           if (((x.id+"\b").replace("_c\b", "_e") not in all_essential_metabolites)) and
                              (((x.id+"\b").replace("_p\b", "_e") not in all_essential_metabolites))]
all_input_metabolite_ids = []
all_output_metabolite_ids = []
all_inout_metabolite_ids_mapping: Dict[str, str] = {}
for periplasmic_metabolite in periplasmic_metabolites:
    metabolite_id = periplasmic_metabolite.id
    all_input_metabolite_ids.append(metabolite_id)
    all_output_metabolite_ids.append(metabolite_id)
    all_inout_metabolite_ids_mapping[metabolite_id] = (metabolite_id+"\b").replace("_p\b", "").replace("_c\b", "")

combined_input_metabolite_ids = list(set(all_input_metabolite_ids + essential_in_metabolites + essential_out_metabolites))
combined_output_metabolite_ids = list(set(all_output_metabolite_ids + essential_out_metabolites + essential_in_metabolites))
combined_metabolite_ids_mapping = {**all_inout_metabolite_ids_mapping, **essential_metabolite_mapping}

# Define pyredcom SingleModel
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
potential_product_metabolite_ids = [(x.id+"\b").replace("_p\b", "").replace("_c\b", "") for x in potential_product_metabolites]

community = Community(
    single_models=[ecoli_1, ecoli_2],
    exchange_compartment_id="exchg",
    exchange_reaction_id_prefix="EX_C_",
    input_metabolite_ids=[(x+"\b").replace("_e\b", "") for x in essential_in_metabolites],
    output_metabolite_ids=[(x+"\b").replace("_e\b", "").replace("\b", "") for x in essential_out_metabolites+potential_product_metabolite_ids]
)
community_model = generate_community_cobra_model(community)


dG0_data_dict = json_load("./pyredcom/publication_run/dG0_data_JSON/dG0.json")


thermodynamically_excluded_metabolites = ["h2o", "h", "pi"]
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
            if reaction.id.endswith("_to_"+(essential_metabolite+"\b").replace("_e\b", "")):
                is_essential = True
                break
        if is_essential:
            continue
        reaction.lower_bound = 0
        reaction.upper_bound = 0
    elif reaction.id.startswith("EX_C"):
        is_essential = False
        for essential_metabolite in all_essential_metabolites:
            if reaction.id.endswith((essential_metabolite+"\b").replace("_e\b", "_exchg")):
                is_essential = True
                break
        if is_essential:
            continue
        reaction.lower_bound = 0
        reaction.upper_bound = 0

print("Done!")
print("===\nRedCom FBA with community model:")
with community_model:
    redcom_fba(community_model, .1)

# Store model as SBML :D
cobra.io.write_sbml_model(community_model, "./pyredcom/publication_run/pyredcom_publication_sbmls_and_dG0_jsons/ecgsDouble__dG0_ARB_no_h2o_no_pi__metabolites_periplasmic_cytosolic.xml")

json_write("./pyredcom/publication_run/pyredcom_publication_sbmls_and_dG0_jsons/ecgsDouble__dG0_ARB_no_h2o_no_pi__metabolites_periplasmic_cytosolic.json", dG0_data_dict)


def run():
    """Function used for calling with "publication_run.py"""
    pass
