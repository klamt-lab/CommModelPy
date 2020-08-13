# This scripts calls all runs and submodules which were used for the generation of
# the community models used in commodepy's publication.
#
# The strange format using "run" is used in order to let the scripts run within
# the commodepy package since these scripts use commodepy for themselves.


from commodepy.publication_run.dG0_data_conversion.dG0_format_to_JSON import run as run_part_1
from commodepy.publication_run.ecgs_load_and_save_with_cobrapy import run as run_part_2
from commodepy.publication_run.ecolicore2_load_and_save_with_cobrapy import run as run_part_3
from commodepy.publication_run.ecgsDouble__dG0_ARB_no_h2o_no_pi__exchangeable_metabolites_periplasmic_cytosolic import run as run_part_4
from commodepy.publication_run.ecolicore2double__dG0_ARB_no_h2o_no_pi__exchangeable_metabolites_periplasmic_cytosolic import run as run_part_5
from commodepy.publication_run.ecolicore2triple__dG0_ARB_no_h2o_no_pi__exchangeable_metabolites_periplasmic_cytosolic import run as run_part_6
from commodepy.publication_run.iML1515double__dG0_ARB_no_h2o_no_pi__exchangeable_metabolites_periplasmic_cytosolic import run as run_part_7

run_part_1()
run_part_2()
run_part_3()
run_part_4()
run_part_5()
run_part_6()
run_part_7()
