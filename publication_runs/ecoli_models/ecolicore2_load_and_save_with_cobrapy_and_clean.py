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
"""This scripts loads and cleans the EcoliCore2 SBML as given in its publication (Hädicke & Klamt, 2017) and saves it again with cobrapy.

This is done in order to gain a version which is not altered by cobrapy while loading it.
The "cleaning" contains steps such as setting 1000 flux bounds to inf.

References:
Hädicke, O., & Klamt, S. (2017). EColiCore2: a reference network model of the central metabolism
of Escherichia coli and relationships to its genome-scale parent model.
Scientific reports, 7, 39647.
"""
import cobra

print("=>Loading and saving of EcoliCore2 while cleaning up some wrong reaction and metabolite ID parts")

print("Loading original EcoliCore2 SBML as given in its publication...")
model = cobra.io.read_sbml_model("publication_runs/ecoli_models/original_sbml_models/ecolicore2.xml")

print("Set -1000/1000 bounds to -inf/inf...")
for reaction in model.reactions:
    if reaction.upper_bound >= 1000:
        reaction.upper_bound = float("inf")
    if reaction.lower_bound <= -1000:
        reaction.lower_bound = -float("inf")

print("Cleaning reaction ID parts...")
for reaction in model.reactions:
    if "_DASH_" in reaction.id:
        reaction.id = reaction.id.replace("_DASH_", "__")
    if "_LPAREN_" in reaction.id:
        reaction.id = reaction.id.replace("_LPAREN_", "_")
    if "_RPAREN_" in reaction.id:
        reaction.id = reaction.id.replace("_RPAREN_", "")

print("Delete boundary-condition-free exchange metabolites...")
metabolite_ids = [x.id for x in model.metabolites]
for metabolite_id in metabolite_ids:
    if metabolite_id.endswith("_ex"):
        metabolite_instance = model.metabolites.get_by_id(metabolite_id)
        model.remove_metabolites([metabolite_instance])

print("Remove exchange metabolite reactions...")
reaction_ids = [x.id for x in model.reactions]
for reaction_id in reaction_ids:
    reaction = model.reactions.get_by_id(reaction_id)
    if reaction.metabolites == {}:
        model.remove_reactions([reaction])
        continue

print("Rename wrong metabolite name parts...")
for metabolite in model.metabolites:
    if "_DASH_" in metabolite.id:
        metabolite.id = metabolite.id.replace("_DASH_", "__")
    if "_LPAREN_" in metabolite.id:
        metabolite.id = metabolite.id.replace("_LPAREN_", "_")
    if "_RPAREN_" in metabolite.id:
        metabolite.id = metabolite.id.replace("_RPAREN_", "")

print("Deactivate all C sources except of D-glucose...")
model.reactions.EX_succ_e.lower_bound = 0
model.reactions.EX_glyc_e.lower_bound = 0
model.reactions.EX_ac_e.lower_bound = 0

print("Test cleaned model with FBA...")
with model:
    print("Single model FBA solution:")
    fba_solution = model.optimize()
    print(model.summary())
    for reaction in model.reactions:
        if not reaction.id.startswith("EX_"):
            continue
        if fba_solution.fluxes[reaction.id] != 0:
            print(f"{reaction.id}: {fba_solution.fluxes[reaction.id] }")
    print("~~~")

print("Delete redundant biomass metabolite...")
model.remove_metabolites([model.metabolites.get_by_id("Biomass")])

print("Saving SBML of cleaned EcoliCore2 model...")
cobra.io.write_sbml_model(model, "./publication_runs/ecoli_models/original_sbml_models_in_cleaned_form/ecolicore2_loaded_and_saved_by_cobrapy_cleaned.xml")

print("Done!")
print("")

for reaction in model.reactions:
    if (reaction.lower_bound < 0) and (reaction.id.startswith("EX_")):
        print(reaction.id)
