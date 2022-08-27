#!/usr/bin/python
import sys
#It is important to set the directory, where  the subdirectory python_utilities is located 
# in the python  path thi is done below
sys.path.append("..")
import os
from python_utilities import GRASTestor
import numpy as np
import pylab as pl
from python_utilities import hepunit as unit
from python_utilities.material import MaterialManager
from python_utilities.GRASAlongStepAnalysisManager import GRASAlongStepAnalysisManager

#Create the testor instance
###########################
theTestor=GRASTestor.GRASTestor()
theMatManager=MaterialManager()
theAnalysisManager=GRASAlongStepAnalysisManager()


def run_test(name_test,target_material,thickness_mat=0.01*unit.mm,
                  phys_list_components=["rmc_em_standard"],
                  prim_particle="e-",emin=unit.keV,emax=10.*unit.MeV,nb_evt=1e4,
                  random_seed=True,with_meshing=False,print_modulo=None,
                  E0=.5*unit.MeV, alpha=1.,type_spectrum="POWER",
                  cut_in_range=0.05*unit.mm):
    theTestor.add_along_step_analysis_test(name_test,thickness_mat,target_material,
                  phys_list_components, prim_particle, emin, emax, nb_evt,
                  random_seed, with_meshing, print_modulo,
                  E0, alpha, type_spectrum, cut_in_range)
    theTestor.perform_a_test(name_test,name_test)


def analysis_electron_spectra(test_name,label,emin,emax,nbins,eprim_min,eprim_max,
                                l_min,l_max,is_fwd=True):
    loc_dir=os.getcwd()
    file_dir=os.path.dirname(__file__)
    if (file_dir != ""):
        os.chdir(file_dir)
    os.chdir("../test_results")
    i=0
    root_file="%s/%s.root" %(test_name,test_name)
    f_vec,emin_vec,emax_vec=theAnalysisManager.GetElectronAlongStepEkin(root_file,nbins,emin,emax,
                                eprim_min,eprim_max,l_min,l_max,is_forward=is_fwd)
        
    evec=0.5*(emin_vec+emax_vec)
    pl.semilogy(evec,f_vec/np.sum(f_vec),label=label)
       
                    
        
emin=8.5
emax=11.5       
material="Tantalum"
thickness_mat=0.16*unit.mm
"""
run_test("AlongForward_8.5_11.5MeV_%s_%.4emm_cut0.1" %(material,thickness_mat/unit.mm),material,thickness_mat,
                  phys_list_components=["rmc_em_standard"],
                  prim_particle="e-",emin=emin,emax=emax,nb_evt=1e6,
                  random_seed=True,with_meshing=False,print_modulo=None,
                  alpha=-1.,type_spectrum="POWER",
                  cut_in_range=0.1*unit.mm)

run_test("AlongAdjoint_8.5_11.5MeV_%s_%.4emm_cut0.1" %(material,thickness_mat/unit.mm),material,thickness_mat,
                  phys_list_components=["rmc_em_standard"],
                  prim_particle="adj_e-",emin=emin,emax=emax,nb_evt=1e6,
                  random_seed=True,with_meshing=False,print_modulo=None,
                  alpha=-1.,type_spectrum="POWER",
                  cut_in_range=0.1*unit.mm)

"""

emin=9.
emax=10.
eprim=10.
n=20
lmin=0.19
lmax=0.4
analysis_electron_spectra("AlongAdjoint_8.5_11.5MeV_Tantalum_1.6000e-01mm_cut0.1","Adj",emin,emax,n,eprim*0.99,eprim,
                                lmin,lmax,is_fwd=False)

analysis_electron_spectra("AlongForward_8.5_11.5MeV_Tantalum_1.6000e-01mm_cut0.1","Fwd",emin,emax,n,eprim*0.99,eprim,
                                lmin,lmax,is_fwd=True)



pl.legend()
pl.show()

