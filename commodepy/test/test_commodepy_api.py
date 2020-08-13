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
"""pytest test range for commodepy's functionalities."""
# IMPORT SECTION
# External modules
import cobra
import cobra.test
import copy
# Internal modules, set so that pytest can work with them
from .. import commodepy


def test_commodepy_api():
    """Fails if anything goes wrong xD. Checks if there were any breaking API changes."""
    #  Let us create a community of two textbook E. coli instances :D
    textbook_model = cobra.test.create_test_model("textbook")
    textbook_model.optimize()
    print(textbook_model.summary())

    single_model_1 = commodepy.SingleModel(
        cobra_model=textbook_model,
        species_abbreviation="ecoli1",
        objective_reaction_id="Biomass_Ecoli_core",
        exchange_reaction_id_prefix="EX_",
        input_metabolite_ids=["glc__D_e", "nh4_e", "o2_e", "pi_e", "h2o_e", "h_e"],
        output_metabolite_ids=["ac_e", "co2_e", "h_e", "h2o_e"],
        model_metabolite_to_exchange_id_mapping={
            "glc__D_e": "glc__D",
            "nh4_e": "nh4",
            "o2_e": "o2",
            "pi_e": "pi",
            "h2o_e": "h2o",
            "h_e": "h",
            "ac_e": "ac",
            "co2_e": "co2",
        }
    )

    single_model_2 = copy.deepcopy(single_model_1)
    single_model_2.species_abbreviation = "ecoli2"

    community = commodepy.Community(
        single_models=[single_model_1, single_model_2],
        exchange_compartment_id="exchg",
        exchange_reaction_id_prefix="EX_C_",
        input_metabolite_ids=["glc__D", "nh4", "o2", "pi", "h2o", "h"],
        output_metabolite_ids=["ac", "co2", "h", "h2o"]
    )

    community_model = commodepy.generate_community_cobra_model(community)
    print("")
    commodepy.redcom_fba(community_model, .1)
    print("")
    commodepy.minimal_species_redcom_fba(community_model, .1)
