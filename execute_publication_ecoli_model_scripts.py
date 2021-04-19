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
"""This scripts calls all runs and submodules which were used for the generation of the community models used in commmodelpy's publication.

The strangely looking format using import only is used in order to let the scripts run within
the commmodelpy package since these scripts use commmodelpy for themselves.

More details on every step of this run can be found in the respective imported scripts.
"""

import publication_runs.ecoli_models.iML1515doubleTC_generation
import publication_runs.ecoli_models.ecolicore2_load_and_save_with_cobrapy_and_clean
import publication_runs.ecoli_models.iML1515_load_and_save_with_cobrapy_and_clean
import publication_runs.ecoli_models.create_bigg_id_to_metanetx_mapping_json_from_iML
import publication_runs.ecoli_models.create_reaction_id_to_dG0_mapping_json_for_ecolicore2
import publication_runs.ecoli_models.create_reaction_id_to_dG0_mapping_json_for_iML1515
import publication_runs.ecoli_models.ecolicore2double_generation
import publication_runs.ecoli_models.ecolicore2triple_generation
import publication_runs.ecoli_models.iML1515double_generation
import publication_runs.ecoli_models.print_dG0_statistics_for_publication_models

import publication_runs.ecoli_models.ecolicore2comp_load_and_save_with_cobrapy_and_clean
import publication_runs.ecoli_models.ecolicore2compDouble_generation
import publication_runs.ecoli_models.ecolicore2compDoubleTC_generation
import publication_runs.ecoli_models.ecolicore2doubleTC_generation
