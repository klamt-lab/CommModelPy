#!/usr/bin/env python3
#
# Copyright 2021 PSB
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
from commmodelpy.submodules.helper_general import json_write, json_load

model = cobra.io.read_sbml_model("./publication_runs/ecoli_models/publication_sbmls_and_dG0_jsons/ecolicore2double_model.xml")

exchg_reaction_ids = [x.id for x in model.reactions if x.id.startswith("EXCHG_")]
for exchg_reaction_id in exchg_reaction_ids:
    reaction = model.reactions.get_by_id(exchg_reaction_id)
    if ("_adp_" in reaction.id) or ("_pi_" in reaction.id) or ("_h_" in reaction.id) or ("_h2o_" in reaction.id):
        continue
    if "ecoli1" in reaction.id:
        id_addition = "ecoli1"
    else:
        id_addition = "ecoli2"
    if ((reaction.lower_bound == 0) and (reaction.upper_bound == 0)):
        new_pseudo_metabolite = cobra.Metabolite(id=reaction.id+"_PSEUDOMET_REV", compartment="exchg")
        reaction.add_metabolites({new_pseudo_metabolite: -1})
        new_pseudo_reaction_1 = cobra.Reaction(id=reaction.id.replace("EXCHG_", "").replace("ecoli1_", "").replace("ecoli2_", "")+"_PSEUDOREAC1_REV_"+id_addition,
        lower_bound=0,
        upper_bound=float("inf"))
        new_pseudo_reaction_1.add_metabolites({
            new_pseudo_metabolite: 1,
        })
        new_pseudo_reaction_2 = cobra.Reaction(id=reaction.id.replace("EXCHG_", "").replace("ecoli1_", "").replace("ecoli2_", "")+"_PSEUDOREAC2_REV_"+id_addition,
        lower_bound=0,
        upper_bound=float("inf"))
        new_pseudo_reaction_2.add_metabolites({
            new_pseudo_metabolite: -1,
            model.metabolites.get_by_id("h2o_c_"+id_addition): -1/3,
            model.metabolites.get_by_id("atp_c_"+id_addition): -1/3,
            model.metabolites.get_by_id("adp_c_"+id_addition): 1/3,
            model.metabolites.get_by_id("h_c_"+id_addition): 1/3,
            model.metabolites.get_by_id("pi_c_"+id_addition): 1/3,
        })
        model.add_reactions([new_pseudo_reaction_1, new_pseudo_reaction_2])
    if ((reaction.lower_bound == 0) and (reaction.upper_bound == 0)):
        new_pseudo_metabolite = cobra.Metabolite(id=reaction.id+"_PSEUDOMET_FWD", compartment="exchg")
        reaction.add_metabolites({new_pseudo_metabolite: 1})
        new_pseudo_reaction_1 = cobra.Reaction(id=reaction.id.replace("EXCHG_", "").replace("ecoli1_", "").replace("ecoli2_", "")+"_PSEUDOREAC1_FWD_"+id_addition,
        lower_bound=0,
        upper_bound=float("inf"))
        new_pseudo_reaction_1.add_metabolites({
            new_pseudo_metabolite: 1
        })
        new_pseudo_reaction_2 = cobra.Reaction(id=reaction.id.replace("EXCHG_", "").replace("ecoli1_", "").replace("ecoli2_", "")+"_PSEUDOREAC2_FWD_"+id_addition,
        lower_bound=0,
        upper_bound=float("inf"))
        new_pseudo_reaction_2.add_metabolites({
            new_pseudo_metabolite: -1,
            model.metabolites.get_by_id("h2o_c_"+id_addition): -1/3,
            model.metabolites.get_by_id("atp_c_"+id_addition): -1/3,
            model.metabolites.get_by_id("adp_c_"+id_addition): 1/3,
            model.metabolites.get_by_id("h_c_"+id_addition): 1/3,
            model.metabolites.get_by_id("pi_c_"+id_addition): 1/3,
        })
        model.add_reactions([new_pseudo_reaction_1, new_pseudo_reaction_2])

print("Test TC model with FBA...")
with model:
    model.reactions.get_by_id("EX_C_glc__D_exchg").lower_bound = -10
    print("TC model FBA solution:")
    fba_solution = model.optimize()
    print(model.summary())
    for reaction in model.reactions:
        if not reaction.id.startswith("EXCHG_"):
            continue
        if fba_solution.fluxes[reaction.id] != 0:
            print(f"{reaction.id}: {fba_solution.fluxes[reaction.id] }")
    print("~~~")

# Store model as SBML
cobra.io.write_sbml_model(
    model, "./publication_runs/ecoli_models/publication_sbmls_and_dG0_jsons/ecolicore2doubleTC_model.xml")

# Write new dG0 data
dG0_data_dict = json_load("./publication_runs/ecoli_models/publication_sbmls_and_dG0_jsons/ecolicore2double_dG0.json")
json_write("./publication_runs/ecoli_models/publication_sbmls_and_dG0_jsons/ecolicore2doubleTC_dG0.json", dG0_data_dict)
