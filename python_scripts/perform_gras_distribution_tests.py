#!/usr/bin/python
import sys
#It is important to set the directory, where  the subdirectory python_utilities is located 
# in the python  path thi is done below
sys.path.append("..")
from python_utilities import GRASTestor




#Create the testor instance
###########################
theTestor=GRASTestor.GRASTestor()
theTestor.add_tests_from_GRAS_distribution()
theTestor.add_examples_from_GRAS_distribution()


#Set the maximum number of test to run ion parallel
#################################################
theTestor.set_max_number_running_tests(3)

#Uncomment/comment the following  to run/not run  all tests
#theTestor.perform_all_tests()


#Run of some specific tests (selection by uncomment/commenting the following lines)
theTestor.perform_a_test("ion_qmd")
#theTestor.perform_a_test("ion_incl_abla")
#theTestor.perform_a_test("ion_shielding")
#theTestor.perform_a_test("analysis_LET_with_interaction")
#theTestor.perform_a_test("example_new_normalisation_none")
#theTestor.perform_a_test("example_new_normalisation_perevent")
#theTestor.perform_a_test("example_new_normalisation_geofactor")
#theTestor.perform_a_test("example_rmc_point_detector")
