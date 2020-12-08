import cobra
import copy
from commodelpy.commodelpy import Community, SingleModel, generate_community_model_with_no_growth
from commodelpy.submodules.helper_general import json_write
from typing import Dict

S = cobra.Metabolite(id="S_c", compartment="c")
A = cobra.Metabolite(id="A_c", compartment="c")
B = cobra.Metabolite(id="B_c", compartment="c")
C = cobra.Metabolite(id="C_c", compartment="c")
P = cobra.Metabolite(id="P_c", compartment="c")
X = cobra.Metabolite(id="X_c", compartment="c")

nil_to_S = cobra.Reaction(
    id="EX_A", lower_bound=0, upper_bound=1000)
S_to_A = cobra.Reaction(id="S_to_A", lower_bound=0, upper_bound=float("inf"))
A_to_B = cobra.Reaction(id="A_to_B", lower_bound=0, upper_bound=float("inf"))
B_to_C = cobra.Reaction(id="B_to_C", lower_bound=0, upper_bound=float("inf"))
C_to_P = cobra.Reaction(id="C_to_P", lower_bound=0, upper_bound=float("inf"))
P__to_nil = cobra.Reaction(
    id="EX_P", lower_bound=0, upper_bound=float("inf"))

nil_to_S.add_metabolites({
    S: 1
})
S_to_A.add_metabolites({
    S: -1,
    A: 1,
    X: 1,
})
A_to_B.add_metabolites({
    A: -1,
    B: 1,
    X: -1,
})
B_to_C.add_metabolites({
    B: -1,
    C: 1,
    X: 1,
})
C_to_P.add_metabolites({
    C: -1,
    P: 1,
    X: -1,
})
P__to_nil.add_metabolites({
    P: -1,
})

single_model = cobra.Model(id_or_model="toy_model")
single_model.add_reactions([
    nil_to_S,
    S_to_A,
    A_to_B,
    B_to_C,
    C_to_P,
    P__to_nil
])
single_model.objective = "EX_P"
with single_model:
    single_model.optimize()
    print(single_model.summary())

model_1 = SingleModel(
    cobra_model=single_model,
    species_abbreviation="species1",
    objective_reaction_id="C_to_P",
    exchange_reaction_id_prefix="EX_",
    input_metabolite_ids=["S_c", "A_c", "B_c", "C_c", "P_c"],
    output_metabolite_ids=["A_c", "B_c", "C_c", "P_c"],
    model_metabolite_to_exchange_id_mapping={
        "S_c": "S",
        "A_c": "A",
        "B_c": "B",
        "C_c": "C",
        "P_c": "P"
    }
)
model_2 = copy.deepcopy(model_1)
model_2.species_abbreviation = "species2"

community = Community(
    single_models=[model_1, model_2],
    exchange_compartment_id="exchg",
    exchange_reaction_id_prefix="EX_C_",
    input_metabolite_ids=["S"],
    output_metabolite_ids=["A", "B", "C", "P"]
)
community_model = generate_community_model_with_no_growth(
    community, {"species1": 0.5, "species2": 0.5})

for reaction in community_model.reactions:
    if reaction.id.startswith("EX_C_"):
        if (not reaction.id.startswith("EX_C_P_")) and (not reaction.id.startswith("EX_C_S_")):
            reaction.lower_bound = 0
            reaction.upper_bound = 0
    if reaction.id.startswith("EXCHG_"):
        if (not reaction.id.startswith("EXCHG_species1_P_")) and (not reaction.id.startswith("EXCHG_species2_P_")) \
           and (not reaction.id.startswith("EXCHG_species1_S_")) and (not reaction.id.startswith("EXCHG_species2_S_")):
            reaction.lower_bound = 0
            reaction.upper_bound = 0

print("A")
community_model.objective = "EX_C_P_exchg"
community_model.reactions.EX_C_P_exchg.upper_bound = 1000
solution = community_model.optimize()
print(community_model.summary())
for flux_key in solution.fluxes.keys():
    if abs(solution.fluxes[flux_key]) > 1e-9:
        print(flux_key, solution.fluxes[flux_key])
print("B")
cobra.io.write_sbml_model(community_model, "./publication_runs/toy_model/toymodelDouble.xml")

pre_dictionary: Dict[str, float] = {
    "S_to_A": 4,
    "A_to_B": -5,
    "B_to_C": -5,
    "C_to_P": 4,
}
dG0_dictionary: Dict[str, Dict[str, float]] = {}
for key in pre_dictionary.keys():
    dG0_dictionary[key+"_species1"] = {}
    dG0_dictionary[key+"_species1"]["dG0"] = pre_dictionary[key]
    dG0_dictionary[key+"_species1"]["uncertainty"] = 0
    dG0_dictionary[key+"_species2"] = {}
    dG0_dictionary[key+"_species2"]["dG0"] = pre_dictionary[key]
    dG0_dictionary[key+"_species2"]["uncertainty"] = 0

for reaction in community_model.reactions:
    if reaction.id.startswith("EXCHG_"):
        dG0_dictionary[reaction.id] = {}
        dG0_dictionary[reaction.id]["dG0"] = 0
        dG0_dictionary[reaction.id]["uncertainty"] = 0

json_write("./publication_runs/toy_model/dG0_toymodelDouble.json", dG0_dictionary)
