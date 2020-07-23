# This scripts calls all runs and submodules which were used for the generation of
# the community models used in pyredcom's publication.
#
# The strange format using "run" is used in order to let the scripts run within
# the pyredcom package since these scripts use pyredcom for themselves.

from pyredcom.publication_run.dG0_data_conversion.dG0_format_to_JSON import run as run_part_1
from pyredcom.publication_run.ecgs_load_and_save_with_cobrapy import run as run_part_2
from pyredcom.publication_run.ecolicore2_load_and_save_with_cobrapy import run as run_part_3
from pyredcom.publication_run.ecgsDouble__dG0_ARB_no_h2o_no_pi__exchangeable_metabolites_periplasmic_cytosolic import run as run_part_4
from pyredcom.publication_run.ecolicore2double__dG0_ARB_no_h2o_no_pi__exchangeable_metabolites_periplasmic_cytosolic import run as run_part_5
from pyredcom.publication_run.ecolicore2triple__dG0_ARB_no_h2o_no_pi__exchangeable_metabolites_periplasmic_cytosolic import run as run_part_6


run_part_1()
run_part_2()
run_part_3()
run_part_4()
run_part_5()
run_part_6()
