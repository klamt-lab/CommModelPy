import cobra
from commodelpy.submodules.helper_general import json_load


def print_statistics(title: str, amount: float, total: float):
    print(title + str(amount) + " (" + str(round(100*(amount/total), 2)) + "%)")


def get_nan_statistics(model_name: str, model_path: str, dG0_json_path: str) -> None:
    dG0s = json_load(dG0_json_path)
    model = cobra.io.read_sbml_model(model_path)

    num_with_dG0 = 0
    num_dG0_free_exchange_reactions = 0
    num_community_exchanges_with_dG0_0 = 0
    num_is_dG0_free_transporter = 0
    num_is_deactivated_nan_reaction = 0
    dG0s_split_keys = []
    for key in dG0s.keys():
        if "_" in key:
            new_key = "_".join(key.split("_")[:-1])
        else:
            new_key = key
        dG0s_split_keys.append(new_key)
    for reaction in model.reactions:
        if "_DASH_" in reaction.id:
            reaction.id = reaction.id.replace("_DASH_", "__")
        if "_LPAREN_" in reaction.id:
            reaction.id = reaction.id.replace("_LPAREN_", "_")
        if "_RPAREN_" in reaction.id:
            reaction.id = reaction.id.replace("_RPAREN_", "")

        if reaction.id.startswith("EX_") or reaction.id.startswith("EX_C_"):
            num_dG0_free_exchange_reactions += 1
        elif reaction.id.startswith("EXCHG_"):
            num_community_exchanges_with_dG0_0 += 1
        elif ("_" in reaction.id) and ("_".join(reaction.id.split("_")[:-1]) in dG0s_split_keys):
            num_with_dG0 += 1
        elif reaction.id in dG0s_split_keys:
            num_with_dG0 += 1
        elif ("tpp_" in reaction.id+"_") or ("t1pp_" in reaction.id+"_") or ("tex_" in reaction.id+"_") or reaction.id.endswith("tpp") or "_".join(reaction.id.split("_")[:-1]).endswith("tex"):
            num_is_dG0_free_transporter += 1
        else:
            num_is_deactivated_nan_reaction += 1

    num_reactions = len(model.reactions)
    print("~"+model_name+"~")
    print("Total number metabolites: "+str(len(model.metabolites)))
    print("Total number reactions: "+str(num_reactions))
    print(num_with_dG0+num_community_exchanges_with_dG0_0+num_dG0_free_exchange_reactions+num_is_dG0_free_transporter+num_is_deactivated_nan_reaction)
    print_statistics("Reactions with dG0: ", num_with_dG0, num_reactions)
    print_statistics("EXCHG reactions with dG0=0: ", num_community_exchanges_with_dG0_0, num_reactions)
    print_statistics("EX/EX_C_ reactions with free dG0: ", num_dG0_free_exchange_reactions, num_reactions)
    print_statistics("Selected transport reactions with free dG0: ", num_is_dG0_free_transporter, num_reactions)
    print_statistics("Deactivated NaN reactions: ", num_is_deactivated_nan_reaction, num_reactions)


get_nan_statistics("ecolicore2triple",
                   "publication_runs/ecoli_models/publication_sbmls_and_dG0_jsons/ecolicore2triple_model.xml",
                   "publication_runs/ecoli_models/publication_sbmls_and_dG0_jsons/ecolicore2triple_dG0.json")
get_nan_statistics("ecolicore2double",
                   "publication_runs/ecoli_models/publication_sbmls_and_dG0_jsons/ecolicore2double_model.xml",
                   "publication_runs/ecoli_models/publication_sbmls_and_dG0_jsons/ecolicore2double_dG0.json")
get_nan_statistics("iML1515double",
                   "publication_runs/ecoli_models/publication_sbmls_and_dG0_jsons/iML1515double_model.xml",
                   "publication_runs/ecoli_models/publication_sbmls_and_dG0_jsons/iML1515double_dG0.json")
get_nan_statistics("EColiCore2",
                   "publication_runs/ecoli_models/original_sbml_models/ecolicore2.xml",
                   "publication_runs/ecoli_models/publication_sbmls_and_dG0_jsons/ecolicore2double_dG0.json")
