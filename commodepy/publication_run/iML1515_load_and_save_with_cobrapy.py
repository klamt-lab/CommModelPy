import cobra


model = cobra.io.read_sbml_model("./commodepy/publication_run/original_sbml_models/iML1515.xml")
for reaction in model.reactions:
    if reaction.upper_bound > 1000:
        reaction.upper_bound = 1000
    if reaction.lower_bound < -1000:
        reaction.lower_bound = -1000
cobra.io.write_sbml_model(model, "./commodepy/publication_run/original_sbml_models_loaded_and_saved_by_cobrapy/iML1515_loaded_and_saved_by_cobrapy.xml")


def run():
    """Function used for calling with "publication_run.py"""
    pass
