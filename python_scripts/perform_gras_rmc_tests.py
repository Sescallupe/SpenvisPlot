#!/usr/bin/python
import sys
#It is important to set the directory, where  the subdirectory python_utilities is located 
# in the python  path thi is done below
sys.path.append("..")
from python_utilities import GRASTestor
import numpy as np
from python_utilities import hepunit as unit
from python_utilities.material import MaterialManager
import scipy as scp
from python_utilities import hepunit as unit


#Create the testor instance
###########################
theTestor=GRASTestor.GRASTestor()
theMatManager=MaterialManager()


def read_spenvis_csv_block(file_name,line_start,line_stop):
    file_obj=open(file_name,'r')
    lines=file_obj.readlines()
    matrix=[]
    for line in lines[line_start:line_stop+1]:
        words=line.split(",")
        vec=[]
        for word in words:
            vec+=[float(word)]
        matrix+=[vec]
    return np.array(matrix)
    
def read_gallileo_spectrum(file_name="../inputs/CIRCULAR_ORBIT_23222.00_km_56.00_deg_trapped_flux.txt"):
    res=read_csv_block(file_name,178,206)
    evec=res[:,0]*unit.MeV
    fint_vec=res[:,1]
    fdiff_vec=res[:,2]
    print fdiff_vec
    return evec,fint_vec,fdiff_vec



def run_test(name_test,layer_thicknesses,materials,
                   particle="e-", emin=unit.keV, emax=100.*unit.MeV, nb_evt_fwd=1e5,
                   nb_evt_adjoint=1e5,
                    random_seed=True,print_modulo=None,
                    E0=20.*unit.MeV, alpha=1.,type_spectrum="EXPO",
                    use_brem=True,use_MS=True,use_Compton=True,use_PhotoElec=True,
                    use_Adjoint_Compton=True, use_Adjoint_PhotoElec=True,
                  use_pionisation=True,cut_in_range=0.1*unit.mm,
                  sensitive_layer_index=1,run_adj=True,run_fwd=True):
            
            spec_type=type_spectrum
            evec=None
            fdiff_vec=None
            if type_spectrum=="JUICE":
                spec_type="USER"
                v0=np.array([1.,2.,3.,5.,7.])*0.1*unit.MeV
                v1=np.array([1.,2.,3.,5.,7.])*unit.MeV
                v2=np.array([1.,2.,3.,5.,7.])*10.*unit.MeV
                v3=np.array([1.,2.,3.,5.,7.,10.])*100.*unit.MeV
                evec1=np.concatenate((v0,v1,v2,v3))
                
                fdiff_vec1=np.array([5.55E+08,2.11e08,1.03E+08,4.66E+07,3.09E+07,
                                    1.21E+07,2.17E+06,7.24E+05,
                                    1.53E+05,5.52E+04,1.87E+04,1.90E+03,3.69E+02,
                                    4.69E+01,1.21E+01,2.87E+00,1.77E-01,3.46E-02,
                                    4.45E-03,1.15E-03,2.76E-04])
                evec=.1 *np.power(emax/.1,np.arange(81)/80.)
                fdiff_vec=np.exp(scp.interp(evec,evec1,np.log(fdiff_vec1)))
                
    
            if run_fwd:
                theTestor.add_1D_rmc_test("%s_fwd_run" %(name_test),layer_thicknesses,materials,
                    is_flat=False, particle=particle, 
                    emin=emin, emax=emax, nb_evt=nb_evt_fwd,
                    random_seed=random_seed,with_meshing=False,print_modulo=print_modulo,
                    E0=E0, alpha=alpha,
                    type_spectrum=spec_type,make_adjoint_sim=False,use_brem=use_brem,
                    use_MS=use_MS,use_Compton=use_Compton,use_Adjoint_Compton=use_Adjoint_Compton,
                    use_PhotoElec=use_PhotoElec,use_Adjoint_PhotoElec=use_Adjoint_PhotoElec,
                  use_pionisation=use_pionisation,cut_in_range=cut_in_range,
                  evec=evec,fdiff_vec=fdiff_vec,sensitive_layer_index=sensitive_layer_index)
                theTestor.perform_a_test("%s_fwd_run" %(name_test),name_test)
            if run_adj:
                theTestor.add_1D_rmc_test("%s_adj_run" %(name_test),layer_thicknesses,materials,
                    is_flat=False, particle=particle, 
                    emin=emin, emax=emax, nb_evt=nb_evt_adjoint,
                    random_seed=random_seed,with_meshing=False,print_modulo=print_modulo,
                    E0=E0, alpha=alpha,
                    type_spectrum=spec_type,make_adjoint_sim=True,use_brem=use_brem,
                    use_MS=use_MS,use_Compton=use_Compton,use_Adjoint_Compton=use_Adjoint_Compton,
                    use_PhotoElec=use_PhotoElec,use_Adjoint_PhotoElec=use_Adjoint_PhotoElec,
                  use_pionisation=use_pionisation,cut_in_range=cut_in_range,
                  evec=evec,fdiff_vec=fdiff_vec,
                  sensitive_layer_index=sensitive_layer_index)
                theTestor.perform_a_test("%s_adj_run" %(name_test),name_test)
            
def run_payload_test(name_test,layer_thicknesses,materials,
                   particle="e-", emin=unit.keV, emax=100.*unit.MeV, nb_evt_fwd=1e5,
                   nb_evt_adjoint=1e5,
                    random_seed=True,print_modulo=None,
                    E0=20.*unit.MeV, alpha=1.,type_spectrum="EXPO",
                    use_brem=True,use_MS=True,use_Compton=True,use_PhotoElec=True,
                    use_Adjoint_Compton=True, use_Adjoint_PhotoElec=True,
                  use_pionisation=True,cut_in_range=0.1*unit.mm,
                  sensitive_volume="component1",component_size=0.5*unit.mm):
            
            spec_type=type_spectrum
            evec=None
            fdiff_vec=None
            if type_spectrum=="JUICE":
                spec_type="USER"
                v0=np.array([1.,2.,3.,5.,7.])*0.1*unit.MeV
                v1=np.array([1.,2.,3.,5.,7.])*unit.MeV
                v2=np.array([1.,2.,3.,5.,7.])*10.*unit.MeV
                v3=np.array([1.,2.,3.,5.,7.,10.])*100.*unit.MeV
                evec1=np.concatenate((v0,v1,v2,v3))
                
                fdiff_vec1=np.array([5.55E+08,2.11e08,1.03E+08,4.66E+07,3.09E+07,
                                    1.21E+07,2.17E+06,7.24E+05,
                                    1.53E+05,5.52E+04,1.87E+04,1.90E+03,3.69E+02,
                                    4.69E+01,1.21E+01,2.87E+00,1.77E-01,3.46E-02,
                                    4.45E-03,1.15E-03,2.76E-04])
                evec=.1 *np.power(emax/.1,np.arange(81)/80.)
                fdiff_vec=np.exp(scp.interp(evec,evec1,np.log(fdiff_vec1)))
                
            
            theTestor.add_payload_rmc_test("%s_fwd_run" %(name_test),layer_thicknesses,materials,
                     particle=particle, 
                    emin=emin, emax=emax, nb_evt=nb_evt_fwd,
                    random_seed=random_seed,with_meshing=False,print_modulo=print_modulo,
                    E0=E0, alpha=alpha,
                    type_spectrum=spec_type,make_adjoint_sim=False,use_brem=use_brem,
                    use_MS=use_MS,use_Compton=use_Compton,use_Adjoint_Compton=use_Adjoint_Compton,
                    use_PhotoElec=use_PhotoElec,use_Adjoint_PhotoElec=use_Adjoint_PhotoElec,
                  use_pionisation=use_pionisation,cut_in_range=cut_in_range,
                  evec=evec,fdiff_vec=fdiff_vec,sensitive_volume=sensitive_volume,
                  component_size=component_size)
            
            theTestor.add_payload_rmc_test("%s_adj_run" %(name_test),layer_thicknesses,materials,
                    particle=particle, 
                    emin=emin, emax=emax, nb_evt=nb_evt_adjoint,
                    random_seed=random_seed,with_meshing=False,print_modulo=print_modulo,
                    E0=E0, alpha=alpha,
                    type_spectrum=spec_type,make_adjoint_sim=True,use_brem=use_brem,
                    use_MS=use_MS,use_Compton=use_Compton,use_Adjoint_Compton=use_Adjoint_Compton,
                    use_PhotoElec=use_PhotoElec,use_Adjoint_PhotoElec=use_Adjoint_PhotoElec,
                  use_pionisation=use_pionisation,cut_in_range=cut_in_range,
                  evec=evec,fdiff_vec=fdiff_vec,
                  sensitive_volume=sensitive_volume,
                  component_size=component_size)
            
            theTestor.perform_a_test("%s_fwd_run" %(name_test),name_test)
            theTestor.perform_a_test("%s_adj_run" %(name_test),name_test)

def run_ngrm_test(name_test,
                   particle="e-", emin=unit.keV, emax=100.*unit.MeV, nb_evt_fwd=1e5,
                   nb_evt_adjoint=1e5,
                    random_seed=True,print_modulo=None,
                    E0=20.*unit.MeV, alpha=1.,type_spectrum="EXPO",
                    use_brem=True,use_MS=True,use_Compton=True,use_PhotoElec=True,
                    use_Adjoint_Compton=True, use_Adjoint_PhotoElec=True,
                  use_pionisation=True,splitting_factor_fwd_gammas=10,cut_in_range=0.1*unit.mm,
                  sensitive_volume="GreatBoxBoard1Component",
                  housing_thickness=1.5*unit.mm,
                                   housing_material="Aluminum",
                                   remove_Ta_plate=False,
                                   remove_back_shield=False,run_adj=True,run_fwd=True,suffix_adj_run="",
                                   suffix_fwd_run=""):
            
            spec_type=type_spectrum
            evec=None
            fdiff_vec=None
            if type_spectrum=="JUICE":
                spec_type="USER"
                v0=np.array([1.,2.,3.,5.,7.])*0.1*unit.MeV
                v1=np.array([1.,2.,3.,5.,7.])*unit.MeV
                v2=np.array([1.,2.,3.,5.,7.])*10.*unit.MeV
                v3=np.array([1.,2.,3.,5.,7.,10.])*100.*unit.MeV
                evec1=np.concatenate((v0,v1,v2,v3))
                
                fdiff_vec1=np.array([5.55E+08,2.11e08,1.03E+08,4.66E+07,3.09E+07,
                                    1.21E+07,2.17E+06,7.24E+05,
                                    1.53E+05,5.52E+04,1.87E+04,1.90E+03,3.69E+02,
                                    4.69E+01,1.21E+01,2.87E+00,1.77E-01,3.46E-02,
                                    4.45E-03,1.15E-03,2.76E-04])
                evec=.1 *np.power(emax/.1,np.arange(81)/80.)
                fdiff_vec=np.exp(scp.interp(evec,evec1,np.log(fdiff_vec1)))
            if type_spectrum=="GALLILEO":
                spec_type="USER"
                emax=6.*unit.MeV
                evec1,fint_vec,fdiff_vec1=read_gallileo_spectrum
                evec=.1 *np.power(emax/.1,np.arange(81)/80.)
                fdiff_vec=np.exp(scp.interp(evec,evec1,np.log(fdiff_vec1)))
                
                
            
            theTestor.add_ngrm_rmc_test("%s_fwd_run%s" %(name_test,suffix_fwd_run),
                     particle=particle, 
                    emin=emin, emax=emax, nb_evt=nb_evt_fwd,
                    random_seed=random_seed,with_meshing=False,print_modulo=print_modulo,
                    E0=E0, alpha=alpha,
                    type_spectrum=spec_type,make_adjoint_sim=False,use_brem=use_brem,
                    use_MS=use_MS,use_Compton=use_Compton,use_Adjoint_Compton=use_Adjoint_Compton,
                    use_PhotoElec=use_PhotoElec,use_Adjoint_PhotoElec=use_Adjoint_PhotoElec,
                  use_pionisation=use_pionisation,splitting_factor_fwd_gammas=splitting_factor_fwd_gammas,
                  cut_in_range=cut_in_range,
                  evec=evec,fdiff_vec=fdiff_vec,sensitive_volume=sensitive_volume,
                  housing_thickness=housing_thickness,
                                   housing_material=housing_material,
                                   remove_Ta_plate=remove_Ta_plate,
                                   remove_back_shield=remove_back_shield)
            
            theTestor.add_ngrm_rmc_test("%s_adj_run%s" %(name_test,suffix_adj_run),
                    particle=particle, 
                    emin=emin, emax=emax, nb_evt=nb_evt_adjoint,
                    random_seed=random_seed,with_meshing=False,print_modulo=print_modulo,
                    E0=E0, alpha=alpha,
                    type_spectrum=spec_type,make_adjoint_sim=True,use_brem=use_brem,
                    use_MS=use_MS,use_Compton=use_Compton,use_Adjoint_Compton=use_Adjoint_Compton,
                    use_PhotoElec=use_PhotoElec,use_Adjoint_PhotoElec=use_Adjoint_PhotoElec,
                  use_pionisation=use_pionisation,splitting_factor_fwd_gammas=splitting_factor_fwd_gammas,
                  cut_in_range=cut_in_range,
                  evec=evec,fdiff_vec=fdiff_vec,
                  sensitive_volume=sensitive_volume,
                  housing_thickness=housing_thickness,
                                   housing_material=housing_material,
                                   remove_Ta_plate=remove_Ta_plate,
                                   remove_back_shield=remove_back_shield)
            
            if run_fwd:
                theTestor.perform_a_test("%s_fwd_run%s" %(name_test,suffix_fwd_run),name_test)
            if run_adj:
                theTestor.perform_a_test("%s_adj_run%s" %(name_test,suffix_adj_run),name_test)

def run_laplace_estec_test(name_test,
                   particle="e-", emin=unit.keV, emax=100.*unit.MeV, nb_evt_fwd=1e5,
                   nb_evt_adjoint=1e5,
                    random_seed=True,print_modulo=None,
                    E0=20.*unit.MeV, alpha=1.,type_spectrum="EXPO",
                    use_brem=True,use_MS=True,use_Compton=True,use_PhotoElec=True,
                    use_Adjoint_Compton=True, use_Adjoint_PhotoElec=True,
                  use_pionisation=True,splitting_factor_fwd_gammas=10,
                  cut_in_range=0.1*unit.mm,
                  sensitive_volume="Target_000",
                  run_adj=True,run_fwd=True,adj_run_nb=1,
                                   fwd_run_nb=1):
            
            spec_type=type_spectrum
            evec=None
            fdiff_vec=None
            if type_spectrum=="JUICE":
                spec_type="USER"
                v0=np.array([1.,2.,3.,5.,7.])*0.1*unit.MeV
                v1=np.array([1.,2.,3.,5.,7.])*unit.MeV
                v2=np.array([1.,2.,3.,5.,7.])*10.*unit.MeV
                v3=np.array([1.,2.,3.,5.,7.,10.])*100.*unit.MeV
                evec1=np.concatenate((v0,v1,v2,v3))
                
                fdiff_vec1=np.array([5.55E+08,2.11e08,1.03E+08,4.66E+07,3.09E+07,
                                    1.21E+07,2.17E+06,7.24E+05,
                                    1.53E+05,5.52E+04,1.87E+04,1.90E+03,3.69E+02,
                                    4.69E+01,1.21E+01,2.87E+00,1.77E-01,3.46E-02,
                                    4.45E-03,1.15E-03,2.76E-04])
                evec=.1 *np.power(emax/.1,np.arange(81)/80.)
                fdiff_vec=np.exp(scp.interp(evec,evec1,np.log(fdiff_vec1)))
            
            if type_spectrum=="GALLILEO":
                spec_type="USER"
                emax=6.*unit.MeV
                evec1,fint_vec,fdiff_vec1=read_gallileo_spectrum
                evec=.1 *np.power(emax/.1,np.arange(81)/80.)
                fdiff_vec=np.exp(scp.interp(evec,evec1,np.log(fdiff_vec1)))
                
            theTestor.add_laplace_estec_rmc_test("fwd_run%i" %(fwd_run_nb),
                     particle=particle, 
                    emin=emin, emax=emax, nb_evt=nb_evt_fwd,
                    random_seed=random_seed,with_meshing=False,print_modulo=print_modulo,
                    E0=E0, alpha=alpha,
                    type_spectrum=spec_type,make_adjoint_sim=False,use_brem=use_brem,
                    use_MS=use_MS,use_Compton=use_Compton,use_Adjoint_Compton=use_Adjoint_Compton,
                    use_PhotoElec=use_PhotoElec,use_Adjoint_PhotoElec=use_Adjoint_PhotoElec,
                  use_pionisation=use_pionisation,splitting_factor_fwd_gammas=splitting_factor_fwd_gammas,
                  cut_in_range=cut_in_range,
                  evec=evec,fdiff_vec=fdiff_vec,sensitive_volume=sensitive_volume)
            
            theTestor.add_laplace_estec_rmc_test("adj_run%i" %(adj_run_nb),
                    emin=emin, emax=emax, nb_evt=nb_evt_adjoint,
                    random_seed=random_seed,with_meshing=False,print_modulo=print_modulo,
                    E0=E0, alpha=alpha,
                    type_spectrum=spec_type,make_adjoint_sim=True,use_brem=use_brem,
                    use_MS=use_MS,use_Compton=use_Compton,use_Adjoint_Compton=use_Adjoint_Compton,
                    use_PhotoElec=use_PhotoElec,use_Adjoint_PhotoElec=use_Adjoint_PhotoElec,
                  use_pionisation=use_pionisation,splitting_factor_fwd_gammas=splitting_factor_fwd_gammas,
                  cut_in_range=cut_in_range,
                  evec=evec,fdiff_vec=fdiff_vec,
                  sensitive_volume=sensitive_volume)
            
            if run_fwd:
                theTestor.perform_a_test("fwd_run%i" %(fwd_run_nb),name_test)
            if run_adj:
                theTestor.perform_a_test("adj_run%i" %(adj_run_nb),name_test)
            
           




#Set the maximum number of test to run in parallel
#################################################
theTestor.set_max_number_running_tests(3)

#run_test("test_1D_5mmAl",np.array([5.,10.])*unit.mm,["Aluminum","Silicon"])
#run_test("test_1D_5mmVacuum",np.array([5.,10.])*unit.mm,["Vacuum","Silicon"])

al_thicknesses=np.array([0.5,1.,5.,10.])*unit.mm

E0_vec=np.array([0.25,0.5,1.,5.,10.,20.])*unit.MeV
emax_vec=np.array([1.,2.,5.,20.,50.,100.])*unit.MeV

E0_vec=np.array([10.,20.])*unit.MeV
emax_vec=np.array([50.,100.])*unit.MeV

materials=["Aluminum","Tantalum"]
"""

al_density=theMatManager.get_material_density("Aluminum")
for al_thickness in al_thicknesses[2:]:
    i=0
    for E0 in E0_vec:
        emax=emax_vec[i]
        i+=1
        for material in materials:
            mat_thickness=al_thickness*al_density/theMatManager.get_material_density(material)
            run_test("test_1D_%.4emmeqAl_%s_%.4eMeVE0_big_cut" %(al_thickness/unit.mm,material,E0/unit.MeV),
                     np.array([mat_thickness/unit.mm,10.])*unit.mm,[material,"Silicon"],
                     E0=E0,emax=emax,cut_in_range=10000.*unit.cm)
            

run_test("test_1D_10mm_Al_E020MeV_cut_big",np.array([10.,10.])*unit.mm,
         ["Aluminum","Silicon"],E0=20.*unit.MeV,emax=100.*unit.MeV,
         cut_in_range=10000.*unit.cm)
"""
al_thicknesses=np.array([5.,10.])*unit.mm
al_thicknesses=np.array([5,10.])*unit.mm
al_thicknesses=np.array([0.05,5.])*unit.mm

E0_vec=np.array([1.,5.,10.,20.])*unit.MeV
emax_vec=np.array([7.,10.,20,50.,100.,1000.])*unit.MeV
#emax_vec=np.array([1.,10.,100.,1000.])*unit.MeV
materials=["Aluminum","Tantalum"]
#materials=["Vacuum"]
#al_thicknesses=np.array([5.])*unit.mm
nb_evt_fwd=1e5
nb_evt_adj=1e5
al_density=theMatManager.get_material_density("Aluminum")
for al_thickness in al_thicknesses[1:]:
    for material in materials[1:]:
        mat_thickness=al_thickness
        if (material != "Vacuum"):
            mat_thickness=al_thickness*al_density/theMatManager.get_material_density(material)
        
        l=0
        """
        for emax in emax_vec[-2:-1]:
            
            run_test("test_1D_%.4emmeqAl_%s_MS_io_and_brem_Juice_emax_%.2f" %(al_thickness/unit.mm,material,emax/unit.MeV),
                     np.array([mat_thickness,10.*unit.mm]),[material,"Silicon"],
                     emax=emax,type_spectrum="JUICE", use_brem=True,
                  use_MS=True,use_Compton=False,use_PhotoElec=False,
                  use_pionisation=False,cut_in_range=0.01*unit.mm, nb_evt_fwd=nb_evt_fwd,
                  nb_evt_adjoint=nb_evt_adj,alpha=-1.,E0=1.,emin=.001*unit.MeV,
                  particle="e-")
            
            run_test("test_1D_%.4emmeqAl_%s_Juice_emax_%.2f" %(al_thickness/unit.mm,material,emax/unit.MeV),
                     np.array([mat_thickness,10.*unit.mm]),[material,"Silicon"],
                     emax=emax,type_spectrum="JUICE", use_brem=True,
                  use_MS=True,use_Compton=True,use_PhotoElec=True,
                  use_pionisation=False,cut_in_range=0.01*unit.mm, nb_evt_fwd=nb_evt_fwd,
                  nb_evt_adjoint=nb_evt_adj,alpha=-1.,E0=1.,emin=.001*unit.MeV,
                  particle="e-")
            
            run_test("test_1D_%.4emmeqAl_%s_noMS_Juice_emax_%.2f" %(al_thickness/unit.mm,material,emax/unit.MeV),
                     np.array([mat_thickness,10.*unit.mm]),[material,"Silicon"],
                     emax=emax,type_spectrum="JUICE", use_brem=True,
                  use_MS=False,use_Compton=True,use_PhotoElec=True,
                  use_pionisation=False,cut_in_range=0.01*unit.mm, nb_evt_fwd=nb_evt_fwd,
                  nb_evt_adjoint=nb_evt_adj,alpha=-1.,E0=1.,emin=.001*unit.MeV,
                  particle="e-")
            
          
        
        
        
        l=0
        for E0 in E0_vec:
            emax=emax_vec[l]
            l+=1
            
            run_test("test_1D_%.4emmeqAl_%s_noMS_e0_%.2f_emax_%.2f_new_along" %(al_thickness/unit.mm,material,
                                                                                  E0/unit.MeV,emax/unit.MeV),
                     np.array([mat_thickness,100.*unit.mm]),[material,"Silicon"],
                     emax=emax,type_spectrum="EXPO", use_brem=True,
                  use_MS=False,use_Compton=True,use_PhotoElec=True,
                  use_pionisation=False,cut_in_range=0.01*unit.mm, nb_evt_fwd=nb_evt_fwd,
                  nb_evt_adjoint=nb_evt_adj,alpha=-1.,E0=E0,emin=.001*unit.MeV,
                  particle="e-")
        
            run_test("test_1D_%.4emmeqAl_%s_noMS_e0_%.2f_emax_%.2f" %(al_thickness/unit.mm,material,
                                                                                  E0/unit.MeV,emax/unit.MeV),
                     np.array([mat_thickness,100.*unit.mm]),[material,"Silicon"],
                     emax=emax,type_spectrum="EXPO", use_brem=True,
                  use_MS=False,use_Compton=True,use_PhotoElec=True,
                  use_pionisation=False,cut_in_range=0.01*unit.mm, nb_evt_fwd=nb_evt_fwd,
                  nb_evt_adjoint=nb_evt_adj,alpha=-1.,E0=E0,emin=.001*unit.MeV,
                  particle="e-")
            run_test("test_1D_%.4emmeqAl_%s_e0_%.2f_emax_%.2f" %(al_thickness/unit.mm,material,
                                                                                  E0/unit.MeV,emax/unit.MeV),
                     np.array([mat_thickness,100.*unit.mm]),[material,"Silicon"],
                     emax=emax,type_spectrum="EXPO", use_brem=True,
                  use_MS=True,use_Compton=True,use_PhotoElec=True,
                  use_pionisation=False,cut_in_range=0.01*unit.mm, nb_evt_fwd=nb_evt_fwd,
                  nb_evt_adjoint=nb_evt_adj,alpha=-1.,E0=E0,emin=.001*unit.MeV,
                  particle="e-")
        
            run_test("test_1D_%.4emmeqAl_%s_noCompton_e0_%.2f_emax_%.2f" %(al_thickness/unit.mm,material,
                                                                                  E0/unit.MeV,emax/unit.MeV),
                     np.array([mat_thickness,100.*unit.mm]),[material,"Silicon"],
                     emax=emax,type_spectrum="EXPO", use_brem=True,
                  use_MS=True,use_Compton=False,use_PhotoElec=True,
                  use_pionisation=False,cut_in_range=0.01*unit.mm, nb_evt_fwd=nb_evt_fwd,
                  nb_evt_adjoint=nb_evt_adj,alpha=-1.,E0=E0,emin=.001*unit.MeV,
                  particle="e-")
        
            run_test("test_1D_%.4emmeqAl_%s_noPEEFFECT_e0_%.2f_emax_%.2f" %(al_thickness/unit.mm,material,
                                                                                  E0/unit.MeV,emax/unit.MeV),
                     np.array([mat_thickness,100.*unit.mm]),[material,"Silicon"],
                     emax=emax,type_spectrum="EXPO", use_brem=True,
                  use_MS=True,use_Compton=True,use_PhotoElec=False,
                  use_pionisation=False,cut_in_range=0.01*unit.mm, nb_evt_fwd=nb_evt_fwd,
                  nb_evt_adjoint=nb_evt_adj,alpha=-1.,E0=E0,emin=.001*unit.MeV,
                  particle="e-")
            
            run_test("test_1D_%.4emmeqAl_%s_MS_io_and_brem_e0_%.2f_emax_%.2f" %(al_thickness/unit.mm,material,
                                                                                  E0/unit.MeV,emax/unit.MeV),
                     np.array([mat_thickness,100.*unit.mm]),[material,"Silicon"],
                     emax=emax,type_spectrum="EXPO", use_brem=True,
                  use_MS=True,use_Compton=False,use_PhotoElec=False,
                  use_pionisation=False,cut_in_range=0.01*unit.mm, nb_evt_fwd=nb_evt_fwd,
                  nb_evt_adjoint=nb_evt_adj,alpha=-1.,E0=E0,emin=.001*unit.MeV,
                  particle="e-")
        """
"""
run_test("test_Ta_2mm_compton_noMS_juice_gamma",
        np.array([2.*unit.mm,5.*unit.mm]),["Tantalum","Silicon"],
        emax=1000.*unit.MeV,type_spectrum="JUICE", use_brem=False,
        use_MS=False,
        use_Compton=True,
        use_PhotoElec=False,
        use_Adjoint_Compton=True,
        use_Adjoint_PhotoElec=True,
        use_pionisation=False,cut_in_range=0.01*unit.mm, nb_evt_fwd=1.e6,
        nb_evt_adjoint=1.e5, alpha=-1.,E0=1.*unit.MeV,emin=.001*unit.MeV,
        particle="gamma",
        sensitive_layer_index=1)   
"""




def run_juice_payload_sim():
    run_payload_test("test_payload_juice_big_s_volume_Ta3mm_io_brem_noMS",
                 np.array([3.*unit.mm]),["Tantalum"],
                 emax=1000.*unit.MeV,type_spectrum="JUICE", use_brem=True,
                 use_MS=False,
                 use_Compton=False,
                 use_PhotoElec=False,
                 use_Adjoint_Compton=True,
                 use_Adjoint_PhotoElec=True,
                 use_pionisation=False,cut_in_range=0.01*unit.mm, nb_evt_fwd=1.e7,
                 nb_evt_adjoint=5.e5, alpha=-1.,E0=1.*unit.MeV,emin=.001*unit.MeV,
                 particle="e-",
            sensitive_volume="component1",
            component_size=5*unit.mm)
    
    run_payload_test("test_payload_juice_big_s_volume_Ta3mm_io_brem",
                 np.array([3.*unit.mm]),["Tantalum"],
                 emax=1000.*unit.MeV,type_spectrum="JUICE", use_brem=True,
                 use_MS=True,
                 use_Compton=False,
                 use_PhotoElec=False,
                 use_Adjoint_Compton=True,
                 use_Adjoint_PhotoElec=True,
                 use_pionisation=False,cut_in_range=0.01*unit.mm, nb_evt_fwd=1.e7,
                 nb_evt_adjoint=5.e5, alpha=-1.,E0=1.*unit.MeV,emin=.001*unit.MeV,
                 particle="e-",
            sensitive_volume="component1",
            component_size=5*unit.mm)
    run_payload_test("test_payload_juice_big_s_volume_Ta3mm_io",
                 np.array([3.*unit.mm]),["Tantalum"],
                 emax=1000.*unit.MeV,type_spectrum="JUICE", use_brem=False,
                 use_MS=True,
                 use_Compton=False,
                 use_PhotoElec=False,
                 use_Adjoint_Compton=True,
                 use_Adjoint_PhotoElec=True,
                 use_pionisation=False,cut_in_range=0.01*unit.mm, nb_evt_fwd=1.e7,
                 nb_evt_adjoint=5.e5, alpha=-1.,E0=1.*unit.MeV,emin=.001*unit.MeV,
                 particle="e-",
            sensitive_volume="component1",
            component_size=5*unit.mm)
    run_payload_test("test_payload_juice_big_s_volume_Ta3mm_io_noMS",
                 np.array([3.*unit.mm]),["Tantalum"],
                 emax=1000.*unit.MeV,type_spectrum="JUICE", use_brem=False,
                 use_MS=False,
                 use_Compton=False,
                 use_PhotoElec=False,
                 use_Adjoint_Compton=True,
                 use_Adjoint_PhotoElec=True,
                 use_pionisation=False,cut_in_range=0.01*unit.mm, nb_evt_fwd=1.e7,
                 nb_evt_adjoint=5.e5, alpha=-1.,E0=1.*unit.MeV,emin=.001*unit.MeV,
                 particle="e-",
            sensitive_volume="component1",
            component_size=5*unit.mm)
#run_juice_payload_sim()    
def run_juice_ngrm_sim(nb_run):
    run_ngrm_test("test_ngrm_juice_2mm_Ta_housing_NoElectronSecondaryButDeltaTest",
                 emax=1000.*unit.MeV,type_spectrum="JUICE", use_brem=True,
                 use_MS=True,
                 use_Compton=True,
                 use_PhotoElec=True,
                 use_Adjoint_Compton=True,
                 use_Adjoint_PhotoElec=False,
                 use_pionisation=False,splitting_factor_fwd_gammas=1,
                 cut_in_range=0.01*unit.mm, nb_evt_fwd=5.e7,
                 nb_evt_adjoint=5.e4, alpha=-1.,E0=0.4*unit.MeV,emin=.001*unit.MeV,
                 particle="e-",housing_thickness=2.*unit.mm,
                                   housing_material="Tantalum",
                                   remove_Ta_plate=True,
                                   remove_back_shield=True,
                                   run_adj=True,run_fwd=False,
                                   suffix_adj_run="5e4evts_splitf1_%i" %(nb_run),
                                   suffix_fwd_run="2")
    
    run_ngrm_test("test_ngrm_juice_2mm_Ta_housing_NoElectronSecondaryButDeltaTest",
                 emax=1000.*unit.MeV,type_spectrum="JUICE", use_brem=True,
                 use_MS=True,
                 use_Compton=True,
                 use_PhotoElec=True,
                 use_Adjoint_Compton=True,
                 use_Adjoint_PhotoElec=False,
                 use_pionisation=False,splitting_factor_fwd_gammas=20,
                 cut_in_range=0.01*unit.mm, nb_evt_fwd=5.e7,
                 nb_evt_adjoint=5.e4, alpha=-1.,E0=0.4*unit.MeV,emin=.001*unit.MeV,
                 particle="e-",housing_thickness=2.*unit.mm,
                                   housing_material="Tantalum",
                                   remove_Ta_plate=True,
                                   remove_back_shield=True,
                                   run_adj=True,run_fwd=False,
                                   suffix_adj_run="5e4evts_splitf20_%i" %(nb_run),suffix_fwd_run="2")
    
     
    
    """
    run_ngrm_test("test_ngrm_juice_10mm_Al_housing_NoElectronSecondaryButDelta",
                 emax=1000.*unit.MeV,type_spectrum="JUICE", use_brem=True,
                 use_MS=True,
                 use_Compton=True,
                 use_PhotoElec=True,
                 use_Adjoint_Compton=True,
                 use_Adjoint_PhotoElec=False,
                 use_pionisation=False,cut_in_range=0.01*unit.mm, nb_evt_fwd=2.e7,
                 nb_evt_adjoint=2.e5, alpha=-1.,E0=0.4*unit.MeV,emin=.001*unit.MeV,
                 particle="e-",housing_thickness=10.*unit.mm,
                                   housing_material="Aluminum",
                                   remove_Ta_plate=True,
                                   remove_back_shield=True)
    """

def run_exp_ngrm_sim():
    run_ngrm_test("test_ngrm_exp_new",
                 emax=7.*unit.MeV,type_spectrum="EXPO", use_brem=True,
                 use_MS=True,
                 use_Compton=True,
                 use_PhotoElec=False,
                 use_Adjoint_Compton=True,
                 use_Adjoint_PhotoElec=False,
                 use_pionisation=False,cut_in_range=0.01*unit.mm, nb_evt_fwd=2.e6,
                 nb_evt_adjoint=5.e4, alpha=-1.,E0=0.25*unit.MeV,emin=.001*unit.MeV,
                 particle="e-")
#run_exp_ngrm_sim()
run_juice_ngrm_sim(7)
def run_juice_laplace_estec_sim():
    for i in range(5):
        run_lapace_estec_test("test_laplace_estec_juice_all_processes_%i" %(i+10),
                 emax=1000.*unit.MeV,type_spectrum="JUICE", 
                 use_brem=True,
                 use_MS=True,
                 use_Compton=True,
                 use_PhotoElec=True,
                 use_Adjoint_Compton=True,
                 use_Adjoint_PhotoElec=True,
                 use_pionisation=False,splitting_factor_fwd_gammas=1,
                 cut_in_range=0.01*unit.mm, nb_evt_fwd=5.e6,
                 nb_evt_adjoint=2.e5, alpha=-1.,E0=0.4*unit.MeV,emin=.001*unit.MeV,
                 particle="e-",
                 suffix_adj_run="5e4evts_splitf1_%i" %(nb_run),
                                   suffix_fwd_run="2")
    
   
#run_juice_laplace_estec_sim()
#run_juice_ngrm_sim()
#run_juice_payload_sim()


def run_juice_sim():
    thicknesses=[10]
    material="Tantalum"
    lab_material="Ta"
    for thickness in thicknesses:
        al_density=theMatManager.get_material_density("Aluminum")
        mat_density=theMatManager.get_material_density(material)
        run_test("test_%s_%.2fmmAlEq_0.01mm_cut_AllNoSecondaryElectronsThirdNoAdjointPhotoElec" %(lab_material,thickness),
                 np.array([thickness*unit.mm*al_density/mat_density,.1*unit.mm,10.*unit.mm]),
                 [material,"Silicon","Vacuum"],
                 emax=1000.*unit.MeV,type_spectrum="JUICE", use_brem=True,
                 use_MS=True,
                 use_Compton=True,
                 use_PhotoElec=True,
                 use_Adjoint_Compton=True,
                 use_Adjoint_PhotoElec=False,
                 use_pionisation=False,cut_in_range=0.01*unit.mm, nb_evt_fwd=5.e6,
                 nb_evt_adjoint=1.e5, alpha=-1.,E0=1.*unit.MeV,emin=.001*unit.MeV,
                 particle="e-",
            sensitive_layer_index=1)
        
    
    """
        run_test("test_%s_%.2fmm_io_noMS_juice" %(lab_material,thickness),
                 np.array([thickness*unit.mm,5.*unit.mm]),[material,"Silicon"],
                 emax=1000.*unit.MeV,type_spectrum="JUICE", use_brem=False,
                 use_MS=False,
                 use_Compton=False,
                 use_PhotoElec=False,
                 use_Adjoint_Compton=False,
                 use_Adjoint_PhotoElec=False,
                 use_pionisation=False,cut_in_range=0.01*unit.mm, nb_evt_fwd=1.e6,
                 nb_evt_adjoint=1.e5, alpha=-1.,E0=1.*unit.MeV,emin=.001*unit.MeV,
                 particle="e-",
                 sensitive_layer_index=1)
    
        run_test("test_%s_%.2fmm_io_brem_juice" %(lab_material,thickness),
                 np.array([thickness*unit.mm,5.*unit.mm]),[material,"Silicon"],
                 emax=1000.*unit.MeV,type_spectrum="JUICE", use_brem=True,
                 use_MS=True,
                 use_Compton=False,
                 use_PhotoElec=False,
                 use_Adjoint_Compton=False,
                 use_Adjoint_PhotoElec=False,
                 use_pionisation=False,cut_in_range=0.01*unit.mm, nb_evt_fwd=1.e6,
                 nb_evt_adjoint=1.e5, alpha=-1.,E0=1.*unit.MeV,emin=.001*unit.MeV,
                 particle="e-",
                 sensitive_layer_index=1)
    
    
        run_test("test_%s_%.2fmm_io_brem_noMS_juice" %(lab_material,thickness),
                 np.array([thickness*unit.mm,5.*unit.mm]),[material,"Silicon"],
                 emax=1000.*unit.MeV,type_spectrum="JUICE", use_brem=True,
                 use_MS=False,
                 use_Compton=False,
            use_PhotoElec=False,
            use_Adjoint_Compton=False,
            use_Adjoint_PhotoElec=False,
            use_pionisation=False,cut_in_range=0.01*unit.mm, nb_evt_fwd=1.e6,
            nb_evt_adjoint=1.e5, alpha=-1.,E0=1.*unit.MeV,emin=.001*unit.MeV,
            particle="e-",
            sensitive_layer_index=1) 
     
     
        run_test("test_%s_%.2fmm_noMS_juice" %(lab_material,thickness),
                 np.array([thickness*unit.mm,5.*unit.mm]),["Aluminum","Silicon"],
                 emax=1000.*unit.MeV,type_spectrum="JUICE", use_brem=True,
                 use_MS=False,
                 use_Compton=True,
                 use_PhotoElec=True,
                 use_Adjoint_Compton=True,
                 use_Adjoint_PhotoElec=True,
                 use_pionisation=False,cut_in_range=0.01*unit.mm, nb_evt_fwd=1.e6,
                 nb_evt_adjoint=1.e5, alpha=-1.,E0=1.*unit.MeV,emin=.001*unit.MeV,
                 particle="e-",
                 sensitive_layer_index=1)
     
        run_test("test_%s_%.2fmm_juice" %(lab_material,thickness),
                 np.array([thickness*unit.mm,5.*unit.mm]),["Aluminum","Silicon"],
                 emax=1000.*unit.MeV,type_spectrum="JUICE", use_brem=True,
                 use_MS=True,
                 use_Compton=True,
                 use_PhotoElec=True,
                 use_Adjoint_Compton=True,
                 use_Adjoint_PhotoElec=True,
                 use_pionisation=False,cut_in_range=0.01*unit.mm, nb_evt_fwd=1.e6,
                 nb_evt_adjoint=1.e5, alpha=-1.,E0=1.*unit.MeV,emin=.001*unit.MeV,
                 particle="e-",
                 sensitive_layer_index=1)
     
        run_test("test_%s_%.2fmm_noMS_noAdjointComptonPhoto_juice" %(lab_material,thickness),
                 np.array([thickness*unit.mm,5.*unit.mm]),["Aluminum","Silicon"],
                 emax=1000.*unit.MeV,type_spectrum="JUICE", use_brem=True,
                 use_MS=False,
                 use_Compton=True,
                 use_PhotoElec=True,
                 use_Adjoint_Compton=False,
                 use_Adjoint_PhotoElec=False,
                 use_pionisation=False,cut_in_range=0.01*unit.mm, nb_evt_fwd=1.e6,
                 nb_evt_adjoint=1.e5, alpha=-1.,E0=1.*unit.MeV,emin=.001*unit.MeV,
                 particle="e-",
                 sensitive_layer_index=1) 
     
        run_test("test_%s_%.2fmm_noAdjointComptonPhoto_juice" %(lab_material,thickness),
                 np.array([thickness*unit.mm,5.*unit.mm]),["Aluminum","Silicon"],
                 emax=1000.*unit.MeV,type_spectrum="JUICE", use_brem=True,
                 use_MS=True,
                 use_Compton=True,
                 use_PhotoElec=True,
                 use_Adjoint_Compton=False,
                 use_Adjoint_PhotoElec=False,
                 use_pionisation=False,cut_in_range=0.01*unit.mm, nb_evt_fwd=1.e6,
                 nb_evt_adjoint=1.e5, alpha=-1.,E0=1.*unit.MeV,emin=.001*unit.MeV,
                 particle="e-",
                 sensitive_layer_index=1)  
    """
#run_juice_sim()        
         


def run_exp_sims(E0=10.*unit.MeV,Emin=0.01*unit.MeV,Emax=1.*unit.MeV,
                 thicknesses=[0.5],
                 material="Aluminum",lab_material="Al"):
    for thickness in thicknesses:
        """
        sim_title="test_%s_%.2fmm_io_E0%.2eMeV_Emin%.2eMeV_Emax%.2eMeV" %(lab_material,thickness,E0/unit.MeV,
                                                                              Emin/unit.MeV,Emax/unit.MeV)
        
        run_test(sim_title,
                 np.array([thickness*unit.mm,5.*unit.mm]),[material,"Silicon"],emin=Emin,E0=E0,
                 emax=Emax,type_spectrum="EXPO", use_brem=False,
                 use_MS=True,
                 use_Compton=False,
                 use_PhotoElec=False,
                 use_Adjoint_Compton=False,
                 use_Adjoint_PhotoElec=False,
                 use_pionisation=False,cut_in_range=0.01*unit.mm, nb_evt_fwd=1.e6,
                 nb_evt_adjoint=1.e5, alpha=-1.,
                 particle="e-",
            sensitive_layer_index=1)
        
        """
        sim_title="test_%s_%.2fmm_io_brem_E0%.2eMeV_Emin%.2eMeV_Emax%.2eMeV" %(lab_material,thickness,E0/unit.MeV,
                                                                              Emin/unit.MeV,Emax/unit.MeV)
        run_test(sim_title,
                 np.array([thickness*unit.mm,5.*unit.mm]),[material,"Silicon"],emin=Emin,E0=E0,
                 emax=Emax,type_spectrum="EXPO", use_brem=True,
                 use_MS=True,
                 use_Compton=False,
                 use_PhotoElec=False,
                 use_Adjoint_Compton=False,
                 use_Adjoint_PhotoElec=False,
                 use_pionisation=False,cut_in_range=0.01*unit.mm, nb_evt_fwd=1.e6,
                 nb_evt_adjoint=1.e5, alpha=-1.,
                 particle="e-",
            sensitive_layer_index=1)
        
        """
        sim_title="test_%s_%.2fmm_E0%.2eMeV_Emin%.2eMeV_Emax%.2eMeV" %(lab_material,thickness,E0/unit.MeV,
                                                                              Emin/unit.MeV,Emax/unit.MeV)
        run_test(sim_title,
                 np.array([thickness*unit.mm,5.*unit.mm]),[material,"Silicon"],emin=Emin,E0=E0,
                 emax=Emax,type_spectrum="EXPO", use_brem=True,
                 use_MS=True,
                 use_Compton=True,
                 use_PhotoElec=True,
                 use_Adjoint_Compton=True,
                 use_Adjoint_PhotoElec=True,
                 use_pionisation=False,cut_in_range=0.01*unit.mm, nb_evt_fwd=1.e6,
                 nb_evt_adjoint=1.e5, alpha=-1.,
                 particle="e-",
            sensitive_layer_index=1)
        
        sim_title="test_%s_%.2fmm_io_noMS_E0%.2eMeV_Emin%.2eMeV_Emax%.2eMeV" %(lab_material,thickness,E0/unit.MeV,
                                                                              Emin/unit.MeV,Emax/unit.MeV)
        run_test(sim_title,
                 np.array([thickness*unit.mm,5.*unit.mm]),[material,"Silicon"],emin=Emin,E0=E0,
                 emax=Emax,type_spectrum="EXPO", use_brem=False,
                 use_MS=False,
                 use_Compton=False,
                 use_PhotoElec=False,
                 use_Adjoint_Compton=False,
                 use_Adjoint_PhotoElec=False,
                 use_pionisation=False,cut_in_range=0.01*unit.mm, nb_evt_fwd=1.e6,
                 nb_evt_adjoint=1.e5, alpha=-1.,
                 particle="e-",
            sensitive_layer_index=1)
        
        sim_title="test_%s_%.2fmm_io_brem_noMS_E0%.2eMeV_Emin%.2eMeV_Emax%.2eMeV" %(lab_material,thickness,E0/unit.MeV,
                                                                              Emin/unit.MeV,Emax/unit.MeV)
        run_test(sim_title,
                 np.array([thickness*unit.mm,5.*unit.mm]),[material,"Silicon"],emin=Emin,E0=E0,
                 emax=Emax,type_spectrum="EXPO", use_brem=True,
                 use_MS=False,
                 use_Compton=False,
                 use_PhotoElec=False,
                 use_Adjoint_Compton=False,
                 use_Adjoint_PhotoElec=False,
                 use_pionisation=False,cut_in_range=0.01*unit.mm, nb_evt_fwd=1.e6,
                 nb_evt_adjoint=1.e5, alpha=-1.,
                 particle="e-",
            sensitive_layer_index=1)
   
run_exp_sims(E0=1.*unit.MeV,Emin=0.01*unit.MeV,Emax=1.*unit.MeV,
                 thicknesses=[.1],
                 material="Aluminum",lab_material="Al") 

run_exp_sims(E0=5.*unit.MeV,Emin=1.*unit.MeV,Emax=10.*unit.MeV,
                 thicknesses=[5.],
                 material="Aluminum",lab_material="Al")  

run_exp_sims(E0=50.*unit.MeV,Emin=10.*unit.MeV,Emax=100.*unit.MeV,
                 thicknesses=[5.],
                 material="Aluminum",lab_material="Al")   

"""


