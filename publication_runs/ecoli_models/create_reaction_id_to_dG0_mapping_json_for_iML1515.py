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
"""Calls 'create_reaction_id_to_dG0_mapping_json' and applies it on iML1515."""
# IMPORT SECTION
# External modules
import cobra
# Internal modules
from publication_runs.ecoli_models.create_reaction_id_to_dG0_mapping_json import create_reaction_id_to_dG0_mapping_json

# ACTUAL ROUTINE SECTION
print("=>Generate reaction ID<->dG0 mapping JSON for EcoliCore2")
ecoli_model = cobra.io.read_sbml_model(
    "./publication_runs/ecoli_models/original_sbml_models_in_cleaned_form/iML1515_loaded_and_saved_by_cobrapy_cleaned.xml")
create_reaction_id_to_dG0_mapping_json(
    ecoli_model, "./publication_runs/ecoli_models/model_reaction_id_to_dG0_mapping_jsons/reaction_id_to_dG0_mapping_for_iML1515.json")
print("Done!")
print("")
