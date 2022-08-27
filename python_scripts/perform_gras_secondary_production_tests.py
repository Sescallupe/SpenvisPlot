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
from python_utilities.GRASReactionAnalyseManager import GRASReactionAnalysisManager

#Create the testor instance
###########################
theTestor=GRASTestor.GRASTestor()
theMatManager=MaterialManager()
theAnalysisManager=GRASReactionAnalysisManager()


def run_test(name_test,target_material,
                 process_name="eBrem",
                  phys_list_components=["rmc_em_standard"],
                  prim_particle="e-",emin=unit.keV,emax=10.*unit.MeV,nb_evt=1e6,
                  random_seed=True,with_meshing=False,print_modulo=None,
                  E0=.5*unit.MeV, alpha=1.,type_spectrum="POWER",
                  cut_in_range=0.05*unit.mm,print_diff_CS=False):
    theTestor.add_process_secondary_analysis_test(name_test,1.*unit.mm,target_material,process_name,
                  phys_list_components, prim_particle, emin, emax, nb_evt,
                  random_seed, with_meshing, print_modulo,
                  E0, alpha, type_spectrum, cut_in_range,print_diff_CS=print_diff_CS)
    theTestor.perform_a_test(name_test,name_test)


emin=.001*unit.MeV
emax=10.*unit.MeV
emin=1.*unit.MeV
emax=1000.*unit.MeV

emin=.01*unit.MeV
emax=.1*unit.MeV
emin=1.*unit.MeV
emax=10.*unit.MeV
emin=.001*unit.MeV
emax=.01*unit.MeV

emin=.0001*unit.MeV
emax=1000.*unit.MeV
material="Tantalum"                  

"""

run_test("InvPEEffect_0.0001_1000MeV_%s" %(material),material,process_name="Inv_PEEffect",
                  phys_list_components=["rmc_em_standard"],
                  prim_particle="adj_e-",emin=emin,emax=emax,nb_evt=1e6,
                  random_seed=True,with_meshing=False,print_modulo=None,
                  alpha=-1.,type_spectrum="POWER",
                  cut_in_range=0.0001*unit.mm)

run_test("PEEffect_0.0001_1000MeV_%s" %(material),material,process_name="phot",
                  phys_list_components=["em_standard"],
                  prim_particle="gamma",emin=emin,emax=emax,nb_evt=1e6,
                  random_seed=True,with_meshing=False,print_modulo=None,
                  alpha=-1.,type_spectrum="POWER",
                  cut_in_range=0.0001*unit.mm) 


run_test("Compton_0.0001_1000MeV_%s" %(material),material,process_name="compt",
                  phys_list_components=["em_standard"],
                  prim_particle="gamma",emin=emin,emax=emax,nb_evt=1e6,
                  random_seed=True,with_meshing=False,print_modulo=None,
                  alpha=-1.,type_spectrum="POWER",
                  cut_in_range=0.0001*unit.mm)

run_test("InvCompt1_0.0001_1000MeV_%s" %(material),material,process_name="Inv_Compt1",
                  phys_list_components=["rmc_em_standard"],
                  prim_particle="adj_e-",emin=emin,emax=emax,nb_evt=1e6,
                  random_seed=True,with_meshing=False,print_modulo=None,
                  alpha=-1.,type_spectrum="POWER",
                  cut_in_range=0.0001*unit.mm)
"""
run_test("InvEbrem1_0.0001_1000MeV_%s" %(material),material,process_name="Inv_eBrem1",
                  phys_list_components=["rmc_em_standard"],
                  prim_particle="adj_gamma",emin=emin,emax=emax,nb_evt=1e5,
                  random_seed=True,with_meshing=False,print_modulo=None,
                  alpha=-1.,type_spectrum="POWER",
                  cut_in_range=0.0001*unit.mm)  

run_test("InvEbrem_0.0001_1000MeV_%s" %(material),material,process_name="Inv_eBrem",
                  phys_list_components=["rmc_em_standard"],
                  prim_particle="adj_e-",emin=emin,emax=emax,nb_evt=1e6,
                  random_seed=True,with_meshing=False,print_modulo=None,
                  alpha=-1.,type_spectrum="POWER",
                  cut_in_range=0.0001*unit.mm)

run_test("Ebrem_0.0001_1000MeV_%s" %(material),material,process_name="eBrem",
                  phys_list_components=["rmc_em_standard"],
                  prim_particle="e-",emin=emin,emax=emax,nb_evt=1e6,
                  random_seed=True,with_meshing=False,print_modulo=None,
                  alpha=-1.,type_spectrum="POWER",
                  cut_in_range=0.0001*unit.mm)          
"""
run_test("InvEioni_0.0001_1000MeV_%s" %(material),material,process_name="Inv_eIon",
                  phys_list_components=["rmc_em_standard"],
                  prim_particle="adj_e-",emin=emin,emax=emax,nb_evt=1e6,
                  random_seed=True,with_meshing=False,print_modulo=None,
                  alpha=-1.,type_spectrum="POWER",
                  cut_in_range=0.0001*unit.mm)

run_test("Eioni_0.0001_1000MeV_%s" %(material),material,process_name="eIoni",
                  phys_list_components=["rmc_em_standard"],
                  prim_particle="e-",emin=emin,emax=emax,nb_evt=1e6,
                  random_seed=True,with_meshing=False,print_modulo=None,
                  alpha=-1.,type_spectrum="POWER",
                  cut_in_range=0.0001*unit.mm)
              
run_test("InvEioni1_0.0001_1000MeV_%s" %(material),material,process_name="Inv_eIon1",
                  phys_list_components=["rmc_em_standard"],
                  prim_particle="adj_e-",emin=emin,emax=emax,nb_evt=1e6,
                  random_seed=True,with_meshing=False,print_modulo=None,
                  alpha=-1.,type_spectrum="POWER",
                  cut_in_range=0.0001*unit.mm)

"""


