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
"""This scripts contains the function which generates the reaction ID<->dG0 mapping as JSON."""
# IMPORT SECTION
# External modules
import cobra
from equilibrator_api import ComponentContribution, Q_
from typing import List, Dict
# Internal modules
from commodelpy.submodules.helper_general import json_load, json_write


def create_reaction_id_to_dG0_mapping_json(model: cobra.Model, json_path: str) -> None:
    """Creates a reaction ID<->dG0 mapping using the Equilibrator API.

    This function uses the pregenerated BIGG ID to MetaNetX ID mapping, and
    lets the Equilibrator API calculate the dG0 values for each reaction using
    the MetaNetX IDs for the metabolites.

    Arguments:
    * model: cobra.Model ~ The cobrapy model instance for which the dG0 mapping shall be created.
    * json_path: str ~

    Output:
    * No variable but a JSON file with the mapping at json_path
    """
    bigg_id_to_metanetx_id = json_load("./publication_runs/ecoli_models/bigg_id_to_metanetx_mapping_json/bigg_id_to_metanetx_id_mapping_from_iML1515.json")
    reaction_id_dG0_mapping: Dict[str, Dict[str, float]] = {}
    cc = ComponentContribution()
    cc.p_h = Q_(7.0)
    cc.ionic_strength = Q_("0.25M")
    cc.temperature = Q_("298.15K")
    cc.p_mg = Q_(3.0)
    reaction_ids_without_dG0: List[str] = []
    for reaction in model.reactions:
        print("==dG0 GENERATION ATTEMPT FOR REACTION "+reaction.id+"==")
        # Exclude exchange reactions
        if reaction.id.startswith("EX_"):
            print("INFO: Reaction is identified as exchange reaction.")
            print("      No dG0 for this reaction!")
            reaction_ids_without_dG0.append(reaction.id)
            continue
        # Exclude special transport reactions
        if ("tpp_" in reaction.id+"_") or ("t1pp_" in reaction.id+"_") or ("tex_" in reaction.id+"_") or reaction.id.endswith("tex"):
            # tpp: Facilitated transport or (proton) symport
            # t1pp: Facilitated transport or (proton) symport
            # tex: "Via diffusion"
            print("INFO: Reaction is identified as special transport reaction.")
            print("      No dG0 for this reaction!")
            reaction_ids_without_dG0.append(reaction.id)
            continue
        # Exclude demand reactions
        if (reaction.id.startswith("DM_")):
            print("INFO: Reaction is identified as sink (demand) reaction.")
            print("      No dG0 for this reaction!")
            reaction_ids_without_dG0.append(reaction.id)
            continue

        educt_strings: List[str] = []
        product_strings: List[str] = []
        is_metanetx_id_missing = False
        for metabolite in reaction.metabolites.keys():
            base_metabolite_id = "".join(metabolite.id.split("_")[:-1])

            if base_metabolite_id.endswith("BIO"):
                base_metabolite_id = base_metabolite_id[:-len("BIO")]

            if base_metabolite_id not in bigg_id_to_metanetx_id.keys():
                is_metanetx_id_missing = True
                break

            if reaction.metabolites[metabolite] < 0:
                educt_strings.append(str(-reaction.metabolites[metabolite]) + " " + bigg_id_to_metanetx_id[base_metabolite_id])
            else:
                product_strings.append(str(reaction.metabolites[metabolite]) + " " + bigg_id_to_metanetx_id[base_metabolite_id])
        if is_metanetx_id_missing:
            print("INFO: MetanetX ID missing for " + base_metabolite_id)
            print("      No dG0 can be given for this reaction!")
            reaction_ids_without_dG0.append(reaction.id)
            continue
        educt_string = " + ".join(educt_strings)
        product_string = " + ".join(product_strings)
        reaction_string = educt_string + " = " + product_string
        print("Reaction string with MetanetX IDs: " + reaction_string)
        try:
            parsed_reaction = cc.parse_reaction_formula(reaction_string)
        except Exception as e:
            print("INFO: Equilibrator reaction parsing error")
            print(e)
            print("      No dG0 can be given for this reaction!")
            reaction_ids_without_dG0.append(reaction.id)
            continue
        try:
            dG0 = cc.standard_dg_prime(parsed_reaction)
            if dG0.m.s > 10000:
                print("INFO: dG0 standard deviation too high")
                print("      No dG0 can be given for this reaction!")
                reaction_ids_without_dG0.append(reaction.id)
                continue
            else:
                reaction_id_dG0_mapping[reaction.id+"_ecoli1"] = {}
                reaction_id_dG0_mapping[reaction.id+"_ecoli1"]["dG0"] = dG0.m.n
                reaction_id_dG0_mapping[reaction.id+"_ecoli1"]["uncertainty"] = 0
                reaction_id_dG0_mapping[reaction.id+"_ecoli2"] = {}
                reaction_id_dG0_mapping[reaction.id+"_ecoli2"]["dG0"] = dG0.m.n
                reaction_id_dG0_mapping[reaction.id+"_ecoli2"]["uncertainty"] = 0
                reaction_id_dG0_mapping[reaction.id+"_ecoli3"] = {}
                reaction_id_dG0_mapping[reaction.id+"_ecoli3"]["dG0"] = dG0.m.n
                reaction_id_dG0_mapping[reaction.id+"_ecoli3"]["uncertainty"] = 0
            print("dG0 calculation successful \\o/")
            print(f"dG0: {dG0}")
        except Exception as e:
            print("INFO:")
            print(e)
            print("      No dG0 can be given for this reaction!")
            reaction_ids_without_dG0.append(reaction.id)
            continue

        # Membrane-bound reaction corrections according to formula 9 in
        # Hamilton, J. J., Dwivedi, V., & Reed, J. L. (2013).
        # Quantitative assessment of thermodynamic constraints on
        # the solution space of genome-scale metabolic models.
        # Biophysical journal, 105(2), 512-522.
        # Formula is:
        # c_j * F * dPsi - 2.3 * h_j * R * T * dpH
        # c_j = net charge transported from outside to inside
        # h_j = Number of protons transported across membrane
        F = 0.10026  # kJ/mV/mol, or 0.02306; kcal/mV/mol
        dPsi = -130  # mV
        dpH = 0.4  # dimensionless
        R = 8.314e-3  # kJ⋅K⁻1⋅mol⁻1
        T = 298.15  # K
        # c_j and h_j are taken from supplementary table 3
        # of the mentioned publication (Hamilton et al., 2013)
        c_j = 0.0
        h_j = 0.0
        # ATP synthase
        if reaction.id.startswith("ATPS4"):
            c_j = 4.0
            h_j = 4.0
        # NADH dehydrogenases
        elif reaction.id.startswith("NADH16"):
            c_j = -3.5
            h_j = -3.5
        elif reaction.id.startswith("NADH17"):
            c_j = -2.0
            h_j = -2.0
        elif reaction.id.startswith("NADH18"):
            c_j = -2.8
            h_j = -2.8
        # Cytochrome dehydrogenases
        elif reaction.id.startswith("CYTBD"):
            c_j = -2.0
            h_j = -2.0
        elif reaction.id.startswith("CYTBO3"):
            c_j = -2.5
            h_j = -2.5

        if (c_j != 0) and (h_j != 0):
            print("Correcting dG0 for special membrane-bound reaction...")
            dG0_correction = c_j * F * dPsi - 2.3 * h_j * R * T * dpH
            print(f"Correction factor is {dG0_correction}")
            reaction_id_dG0_mapping[reaction.id+"_ecoli1"]["dG0"] += dG0_correction
            reaction_id_dG0_mapping[reaction.id+"_ecoli2"]["dG0"] += dG0_correction
            reaction_id_dG0_mapping[reaction.id+"_ecoli3"]["dG0"] += dG0_correction
            print("New dG0 is " + str(reaction_id_dG0_mapping[reaction.id+"_ecoli1"]["dG0"]))

    num_reactions_without_dG0 = len(reaction_ids_without_dG0)
    print("\n==FINAL STATISTICS==")
    print("No dG0 for the following reactions:")
    print("\n".join(reaction_ids_without_dG0))
    print("Total number: ", str(num_reactions_without_dG0))
    print("Fraction from all reactions: ", num_reactions_without_dG0/len(model.reactions))

    json_write(json_path, reaction_id_dG0_mapping)
