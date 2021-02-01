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
"""
This script generates a BiGG ID to MetaNetX ID mapping in the form
of a JSON. This mapping is later used for the Equilibrator API-based
calculation of dG0 values for reactions in iML1515 and EcoliCore2.
"""
import cobra
from commmodelpy.submodules.helper_general import json_write
from typing import Dict

print("=>Create BIGG ID to MetaNetX mapping using the data given in iML1515")

print("Load iML1515 in used form...")
model = cobra.io.read_sbml_model(
    "./publication_runs/ecoli_models/original_sbml_models_in_cleaned_form/iML1515_loaded_and_saved_by_cobrapy_cleaned.xml")

print("Read out BiGG IDs and associated MetaNetX IDs as given in iML1515's reactions, thereby creating the mapping...")
bigg_id_metanetx_id_mapping: Dict[str, Dict[str, str]] = {}
for metabolite in model.metabolites:
    if ("bigg.metabolite" in metabolite.annotation.keys()) and ("metanetx.chemical" in metabolite.annotation.keys()):
        bigg_id = metabolite.annotation["bigg.metabolite"]
        if bigg_id not in bigg_id_metanetx_id_mapping.keys():
            bigg_id_metanetx_id_mapping[bigg_id] = {}

            if "_" in bigg_id:
                bigg_id_metanetx_id_mapping[bigg_id.replace("_", "")] = {}

        inchi_annotation = metabolite.annotation["metanetx.chemical"]
        if type(inchi_annotation) is list:
            selected_inchi_id = inchi_annotation[0]
        else:
            selected_inchi_id = inchi_annotation

        bigg_id_metanetx_id_mapping[bigg_id] = selected_inchi_id
        if "_" in bigg_id:
            bigg_id_metanetx_id_mapping[bigg_id.replace(
                "_", "")] = selected_inchi_id

print("Save the mapping as JSON...")
json_write("./publication_runs/ecoli_models/bigg_id_to_metanetx_mapping_json/bigg_id_to_metanetx_id_mapping_from_iML1515.json",
           bigg_id_metanetx_id_mapping)

print("Done!")
print("")
