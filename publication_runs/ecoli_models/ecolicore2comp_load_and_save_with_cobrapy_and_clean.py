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
model = cobra.io.read_sbml_model("publication_runs/ecoli_models/original_sbml_models/ecolicore2compressed.xml")

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


print("Delete boundary-condition-free exchange metabolites ending with _ex...")
boundary_free_metabolite_ids = [
    x.id.replace("EX_", "") for x in model.reactions
    if x.id.startswith("EX_") and x.id.endswith("_ex")
]
boundary_free_metabolites = []
for x in boundary_free_metabolite_ids:
    try:
        metabolite = model.metabolites.get_by_id(x)
        boundary_free_metabolites.append(metabolite)
    except KeyError:
        continue
for metabolite in boundary_free_metabolites:
    model.remove_metabolites([metabolite])

print("Rename wrong metabolite name parts...")
for metabolite in model.metabolites:
    if "_DASH_" in metabolite.id:
        metabolite.id = metabolite.id.replace("_DASH_", "__")
    if "_LPAREN_" in metabolite.id:
        metabolite.id = metabolite.id.replace("_LPAREN_", "_")
    if "_RPAREN_" in metabolite.id:
        metabolite.id = metabolite.id.replace("_RPAREN_", "")

print("Delete redundant glucose _ex metabolite...")
model.remove_metabolites([model.metabolites.get_by_id("glc__D_ex")])

print("Delete unused boundary-free metabolite exchange reactions...")
ex_reaction_ids = [x.id for x in model.reactions if x.id.startswith("EX_")]
for ex_reaction_id in ex_reaction_ids:
    if model.reactions.get_by_id(ex_reaction_id).metabolites == {}:
        model.remove_reactions([ex_reaction_id])

print("Add periplasmic metabolite for all metabolties in EX_ reactions ending with _c...")
ex_c_reaction_ids = [x.id for x in model.reactions if x.id.startswith("EX_") and x.id.endswith("_c")]
for ex_c_reaction_id in ex_c_reaction_ids:
    reaction = model.reactions.get_by_id(ex_c_reaction_id)
    metabolite = list(reaction.metabolites.keys())[0]
    new_p_metabolite = cobra.Metabolite(id=metabolite.id.replace("_c", "_p"), compartment="p")
    reaction.add_metabolites({
        new_p_metabolite: 1
    })
    new_ex_reaction = cobra.Reaction(id="EX_"+new_p_metabolite.id,
        lower_bound=reaction.lower_bound,
        upper_bound=reaction.upper_bound)
    new_ex_reaction.add_metabolites ({
        new_p_metabolite: -1
    })
    model.add_reactions([new_ex_reaction])
    reaction.id = "Transport_c_to_p_" + metabolite.id.replace("_c", "")

print("Get all reaction IDs ending with Ex and Up...")
ex_reaction_ids = [x.id for x in model.reactions if x.id.endswith("Ex")]
up_reaction_ids = [x.id for x in model.reactions if x.id.endswith("Up")]
all_exchange_ids = ex_reaction_ids + up_reaction_ids

class ExchangedMetabolite:
    def __init__(self, metabolite_id, ex_reaction_id, up_reaction_id):
        self.metabolite_id = metabolite_id
        self.ex_reaction_id = ex_reaction_id
        self.up_reaction_id = up_reaction_id

exchanged_metabolites = []
for exchange_id in all_exchange_ids:
    ex_reaction_id = ""
    up_reaction_id = ""

    if exchange_id.endswith("Up"):
        up_reaction_id = exchange_id
    elif exchange_id.endswith("Ex"):
        ex_reaction_id = exchange_id
    else:
        print(exchange_id)
        print("Error 1!")
        input()

    reaction = model.reactions.get_by_id(exchange_id)
    reaction_id_start = reaction.id[:2].lower()
    exchanged_metabolite = None
    for metabolite in reaction.metabolites:
        if metabolite.id.startswith(reaction_id_start):
            if exchanged_metabolite != None:
                print("Error 1B!")
                input()
            exchanged_metabolite = metabolite

    current_ids = [x.metabolite_id for x in exchanged_metabolites]
    if exchanged_metabolite.id in current_ids:
        element_index = current_ids.index(exchanged_metabolite.id)
        if ex_reaction_id != "":
            exchanged_metabolites[element_index].ex_reaction_id = ex_reaction_id
        elif up_reaction_id != "":
            exchanged_metabolites[element_index].up_reaction_id = up_reaction_id
        else:
            print("Error 2!")
            input()
    else:
        if ex_reaction_id != "":
            up_reaction_id = ""
        elif up_reaction_id != "":
            ex_reaction_id = ""
        else:
            print("Error 3!")
            input()
        exchanged_metabolites.append(ExchangedMetabolite(
            metabolite_id=exchanged_metabolite.id,
            ex_reaction_id=ex_reaction_id,
            up_reaction_id=up_reaction_id
        ))

print("Create new EX_ metabolites with, if not given, new periplasmic intermediates...")
for exchanged_metabolite in exchanged_metabolites:
    has_ex = exchanged_metabolite.ex_reaction_id != ""
    has_up = exchanged_metabolite.up_reaction_id != ""

    if exchanged_metabolite.metabolite_id.endswith("_p"):
        if has_ex and has_up:
            ex_reaction = model.reactions.get_by_id(exchanged_metabolite.ex_reaction_id)
            up_reaction = model.reactions.get_by_id(exchanged_metabolite.up_reaction_id)

            new_ex_reaction = cobra.Reaction(id="EX_"+exchanged_metabolite.metabolite_id,
                lower_bound=-up_reaction.upper_bound,
                upper_bound=ex_reaction.upper_bound)
            new_ex_reaction.add_metabolites({
                model.metabolites.get_by_id(exchanged_metabolite.metabolite_id): -1
            })
            model.add_reactions([new_ex_reaction])
            model.remove_reactions([
                ex_reaction,
                up_reaction
            ])
        elif has_ex:
            ex_reaction = model.reactions.get_by_id(exchanged_metabolite.ex_reaction_id)
            if len(list(ex_reaction.metabolites.keys())) > 1:
                print("Error A1!")
                input()
            ex_reaction.id = "EX_" + exchanged_metabolite.metabolite_id
        elif has_up:
            up_reaction = model.reactions.get_by_id(exchanged_metabolite.up_reaction_id)
            if len(list(up_reaction.metabolites.keys())) > 1:
                print("Error A2!")
                input()
            up_reaction.id = "EX_" + exchanged_metabolite.metabolite_id

            old_lower_bound = up_reaction.lower_bound
            old_upper_bound = up_reaction.upper_bound
            if exchanged_metabolite.metabolite_id.startswith("glc__"):
                print("A")
            up_reaction.lower_bound = -old_upper_bound
            up_reaction.upper_bound = -old_lower_bound

            up_reaction.add_metabolites({
                model.metabolites.get_by_id(exchanged_metabolite.metabolite_id): -2
            })
        else:
            print("Error Zeta!")
            input()
    elif exchanged_metabolite.metabolite_id.endswith("_c"):
        new_p_metabolite_id = exchanged_metabolite.metabolite_id.replace("_c", "_p")
        new_p_metabolite = cobra.Metabolite(id=new_p_metabolite_id, compartment="p")
        model.add_metabolites(new_p_metabolite)

        new_p_ex_reaction = cobra.Reaction(id="EX_" + new_p_metabolite_id)
        new_p_ex_reaction.add_metabolites({
            new_p_metabolite: -1
        })

        if has_ex and has_up:
            ex_reaction = model.reactions.get_by_id(exchanged_metabolite.ex_reaction_id)
            up_reaction = model.reactions.get_by_id(exchanged_metabolite.up_reaction_id)

            new_p_ex_reaction.lower_bound = -up_reaction.upper_bound
            new_p_ex_reaction.upper_bound = ex_reaction.upper_bound

            ex_reaction.add_metabolites({
                new_p_metabolite: 1
            })
            up_reaction.add_metabolites({
                new_p_metabolite: -1
            })
            ex_reaction.id = "Transport_c_to_p_" + new_p_metabolite.id.replace("_p", "")
            up_reaction.id = "Transport_p_to_c_" + new_p_metabolite.id.replace("_p", "")
        elif has_ex:
            ex_reaction = model.reactions.get_by_id(exchanged_metabolite.ex_reaction_id)

            new_p_ex_reaction.lower_bound = 0
            new_p_ex_reaction.upper_bound = ex_reaction.upper_bound

            ex_reaction.add_metabolites({
                new_p_metabolite: 1
            })
            ex_reaction.id = "Transport_c_to_p_" + new_p_metabolite.id.replace("_p", "")
        elif has_up:
            up_reaction = model.reactions.get_by_id(exchanged_metabolite.up_reaction_id)

            new_p_ex_reaction.lower_bound = -up_reaction.upper_bound

            up_reaction.add_metabolites({
                new_p_metabolite: -1
            })
            up_reaction.id = "Transport_p_to_c_" + new_p_metabolite.id.replace("_p", "")
        else:
            print("Error Beta!")
            input()

        model.add_reactions([new_p_ex_reaction])
    else:
        print("Error Alpha!")
        input()

print("Delete biomass metabolite and associated reaction...")
model.remove_metabolites([model.metabolites.get_by_id("Biomass")])
model.remove_reactions([model.reactions.get_by_id("EX_Biomass")])

print("Deactivate all C sources except of D-glucose...")
model.reactions.EX_succ_p.lower_bound = 0
model.reactions.EX_glyc_p.lower_bound = 0
model.reactions.EX_ac_p.lower_bound = 0
model.reactions.EX_glc__D_p.lower_bound = -10

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

print("Print exchange metabolites...")
in_string = ""
out_string = ""
for reaction in model.reactions:
    if reaction.id.startswith("EX_"):
        string_part = f'"{reaction.id.replace("EX_", "")}",\n'
        if reaction.lower_bound < 0:
            in_string += string_part
        if reaction.upper_bound > 0:
            out_string += string_part
print("In metabolites:")
print(in_string)
print("Out metabolites:")
print(out_string)

print("Saving SBML of cleaned EcoliCore2compressed model...")
cobra.io.write_sbml_model(model, "./publication_runs/ecoli_models/original_sbml_models_in_cleaned_form/ecc2comp_loaded_and_saved_by_cobrapy_cleaned.xml")

print("Done!")
print("")

for reaction in model.reactions:
    if (reaction.lower_bound < 0) and (reaction.id.startswith("EX_")):
        print(reaction.id)
