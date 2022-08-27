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

"""
Select if the check test will plot the csv histograms 
"""
view_plot=True

"""
Select the tests for wich the results will be checked
"""
test_names=["example_new_normalisation_tofluence","analysis_LET_with_interaction"]
test_names=["ion_shielding","ion_qmd"]
test_names=["ion_qmd"]
for name in test_names:
    theTestor.first_look_at_test_results(name,
                                         look_gras_doubles=True,
                                         look_log_file=True,
                                         plot_gras_1Dhistos=view_plot)