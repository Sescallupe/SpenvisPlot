#!/usr/bin/python
import sys
#It is important to set the directory, where  the subdirectory python_utilities is located 
# in the python  path thi is done below
sys.path.append("..")
from python_utilities import GRASTestor
import numpy as np
from python_utilities import hepunit as unit



#Create the testor instance
###########################
theTestor=GRASTestor.GRASTestor()

theTestor.add_1D_test("test_1D",np.array([1.,1.])*unit.mm,["Aluminum","Silicon"],
                    is_flat=False, particle="e-", emin=unit.keV, emax=10.*unit.MeV, nb_evt=1e6,
                    random_seed=True,with_meshing=False,print_modulo=None,
                    E0=.5*unit.MeV, alpha=1.,type_spectrum="EXPO")
theTestor.add_1D_rmc_test("test_1D_rmc_1.8mmTa",np.array([1.8,10.])*unit.mm,["Tantalum","Silicon"],
                    is_flat=False, particle="e-", emin=unit.keV, emax=100.*unit.MeV, nb_evt=1e5,
                    random_seed=True,with_meshing=False,print_modulo=None,
                    E0=20.*unit.MeV, alpha=1.,type_spectrum="EXPO")

theTestor.add_1D_rmc_test("test_1D_rmc_fwd_mode_1.8mmTa",np.array([1.8,10.])*unit.mm,["Tantalum","Silicon"],
                    is_flat=False, particle="e-", emin=unit.keV, emax=100.*unit.MeV, nb_evt=1e5,
                    random_seed=True,with_meshing=False,print_modulo=None,
                    E0=20.*unit.MeV, alpha=1.,type_spectrum="EXPO",make_adjoint_sim=False)


#Set the maximum number of test to run ion parallel
#################################################
theTestor.set_max_number_running_tests(3)

theTestor.perform_a_test("test_1D_rmc_fwd_mode_1.8mmTa","test_1D_rmc_fwd_mode_1.8mmTa")
