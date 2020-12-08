import cobra

community_model_path = "publication_runs/ecoli_models/publication_sbmls_and_dG0_jsons/ecolicore2triple_model.xml"
model = cobra.io.read_sbml_model(community_model_path)

found_exchange_reactions = [
"R_EX_C_r1p_exchg",
"R_EX_C_4pasp_exchg",
"R_EX_C_tyr__L_exchg",
"R_EX_C_h_exchg",
"R_EX_C_nicrnt_exchg",
"R_EX_C_pep_exchg",
"R_EX_C_mal__L_exchg",
"R_EX_C_kdo8p_exchg",
"R_EX_C_phpyr_exchg",
"R_EX_C_skm_exchg",
"R_EX_C_3pg_exchg",
"R_EX_C_chor_exchg",
"R_EX_C_3dhsk_exchg",
"R_EX_C_phe__L_exchg",
"R_EX_C_ppi_exchg",
"R_EX_C_h2o_exchg",
"R_EX_C_indole_exchg",
"R_EX_C_skm5p_exchg",
"R_EX_C_oaa_exchg",
"R_EX_C_pram_exchg",
"R_EX_C_2dda7p_exchg",
"R_EX_C_4adcho_exchg",
"R_EX_C_ala__B_exchg",
"R_EX_C_2pg_exchg",
"R_EX_C_pphn_exchg",
"R_EX_C_fum_exchg",
"R_EX_C_3ig3p_exchg",
"R_EX_C_pran_exchg",
"R_EX_C_kdo_exchg",
"R_EX_C_13dpg_exchg",
"R_EX_C_for_exchg",
"R_EX_C_34hpp_exchg",
"R_EX_C_g3p_exchg",
"R_EX_C_glx_exchg",
"R_EX_C_phthr_exchg",
"R_EX_C_2cpr5p_exchg",
"R_EX_C_f6p_exchg",
"R_EX_C_4per_exchg",
"R_EX_C_trp__L_exchg",
"R_EX_C_2me4p_exchg",
"R_EX_C_h2_exchg",
"R_EX_C_3dhq_exchg",
"R_EX_C_3psme_exchg",
"R_EX_C_dhap_exchg",
"R_EX_C_asp__L_exchg",
"R_EX_C_suchms_exchg"
]

output = ""
for reaction in model.reactions:
    if reaction.id.startswith("EX_C_"):
        if "R_"+reaction.id not in found_exchange_reactions:
            output = output + '"R_' + reaction.id + '"\n'

print(output)
