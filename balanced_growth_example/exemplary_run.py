import cobra
import copy
from commodelpy.commodelpy import Community, SingleModel, create_community_model_with_balanced_growth
from typing import Dict

growth_rate = 0.4
# ecolicore double model
# with internal exchanges for everything in the periplasm
ecoli_model = cobra.io.read_sbml_model(
    "./publication_runs/ecoli_models/original_sbml_models/iML1515.xml")

# Deactivate all C sources except of glucose
ecoli_model.reactions.EX_succ_e.lower_bound = 0
ecoli_model.reactions.EX_glyc_e.lower_bound = 0
ecoli_model.reactions.EX_ac_e.lower_bound = 0

# Test model with FBA
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
all_essential_metabolites = list(
    set(essential_in_metabolites+essential_out_metabolites))
essential_metabolite_mapping: Dict[str, str] = {}
for met in essential_in_metabolites+essential_out_metabolites:
    essential_metabolite_mapping[met] = (met+"\b").replace("_e\b", "")

periplasmic_metabolites = [
    x for x in ecoli_model.metabolites if x.id.endswith("_p")]
periplasmic_metabolites_cut = [(x.id+"\b").replace("_p\b", "")
                               for x in periplasmic_metabolites]
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
    all_inout_metabolite_ids_mapping[metabolite_id] = (
        metabolite_id+"\b").replace("_p\b", "").replace("_c\b", "")

combined_input_metabolite_ids = list(set(
    all_input_metabolite_ids + essential_in_metabolites + essential_out_metabolites))
combined_output_metabolite_ids = list(set(
    all_output_metabolite_ids + essential_out_metabolites + essential_in_metabolites))
combined_metabolite_ids_mapping = {
    **all_inout_metabolite_ids_mapping, **essential_metabolite_mapping}

# Define commodelpy SingleModel
ecoli_1 = SingleModel(
    cobra_model=ecoli_model,
    species_abbreviation="ecoli1",
    objective_reaction_id="BIOMASS_Ec_iML1515_core_75p37M",
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
    input_metabolite_ids=[(x+"\b").replace("_e\b", "")
                          for x in essential_in_metabolites],
    output_metabolite_ids=[(x+"\b").replace("_e\b", "").replace("\b", "")
                           for x in essential_out_metabolites+potential_product_metabolite_ids]
)
community_model = create_community_model_with_balanced_growth(
    community, growth_rate)


cobra.io.write_sbml_model(
    community_model, "./balanced_growth_example/models/community_model.xml")

with community_model:
    community_model.optimize()
    print(community_model.summary())


def run():
    pass
