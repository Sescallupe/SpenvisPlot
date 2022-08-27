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



def analysis_gamma_spectra(test_names,labels,with_print_diff_CS=False):
    loc_dir=os.getcwd()
    file_dir=os.path.dirname(__file__)
    if (file_dir != ""):
        os.chdir(file_dir)
    os.chdir("../test_results")
    i=0
    for test_name in test_names:
        root_file="%s/%s.root" %(test_name,test_name)
        f_vec,emin_vec,emax_vec=theAnalysisManager.GetSecondaryGammaSpectrum(root_file,70,0.01*unit.keV,
                                                     100.*unit.MeV,
                                logE=True,tree_file_name="secondary_tuple_interaction")
        
        evec=0.5*(emin_vec+emax_vec)
        pl.loglog(evec,f_vec/10.,label=labels[i])
        if (with_print_diff_CS):
            log_file_diffCS="%s_diffCS/%s_diffCS.log" %(test_name,test_name)
            fo=open(log_file_diffCS,'r')
            lines=fo.readlines()
            evec=[]
            diff_CS_vec=[]
            for line in lines:
                if "Diff macro CS" in line:
                    words=line.split("Second Kin. energy [MeV]:")[1].split("Diff macro CS . [cm/MeV]:")
                    e=float(words[0])
                    CS=float(words[1])
                    if CS >1.e-5:
                        evec+=[float(words[0])]
                        diff_CS_vec+=[float(words[1])]
            pl.loglog(evec,np.array(diff_CS_vec),label=labels[i])
                    
        
        
        i+=1
    #pl.legend(loc="lower left")
    os.chdir(loc_dir)
def analysis_electron_spectra(test_names,labels,with_print_diff_CS=False,extra_cut=""):
    loc_dir=os.getcwd()
    file_dir=os.path.dirname(__file__)
    if (file_dir != ""):
        os.chdir(file_dir)
    os.chdir("../test_results")
    i=0
    for test_name in test_names:
        root_file="%s/%s.root" %(test_name,test_name)
        f_vec,emin_vec,emax_vec=theAnalysisManager.GetSecondaryElectronSpectrum(root_file,70,0.01*unit.keV,
                                                     100.*unit.MeV,
                                logE=True,tree_file_name="secondary_tuple_interaction",extra_cut=extra_cut)
        
        evec=0.5*(emin_vec+emax_vec)
        pl.semilogx(evec,f_vec/10.,label=labels[i])
        print np.sum((emax_vec-emax_vec)*f_vec)
        if (with_print_diff_CS):
            log_file_diffCS="%s_diffCS/%s_diffCS.log" %(test_name,test_name)
            fo=open(log_file_diffCS,'r')
            lines=fo.readlines()
            evec=[]
            diff_CS_vec=[]
            for line in lines:
                if "Diff macro CS" in line:
                    words=line.split("Second Kin. energy [MeV]:")[1].split("Diff macro CS . [cm/MeV]:")
                    e=float(words[0])
                    CS=float(words[1])
                    if CS >1.e-5:
                        evec+=[float(words[0])]
                        diff_CS_vec+=[float(words[1])]
            pl.loglog(evec,np.array(diff_CS_vec),label=labels[i])
                    
        
        
        i+=1
    #pl.legend(loc="lower left")
    os.chdir(loc_dir)
def analysis_elec_spectra(test_name,label_name,eprim_min,eprim_max,f_expression="1./ekin_prim",nprim=1.e6,extra_cut="",
                                 prim_flux_emin=0.001,prim_flux_emax=10000000.,symbol=""):
    loc_dir=os.getcwd()
    file_dir=os.path.dirname(__file__)
    if (file_dir != ""):
        os.chdir(file_dir)
    os.chdir("../test_results")
    root_file="%s/%s.root" %(test_name,test_name)
    weight_correction="(ekin_prim)*(%s)" %( f_expression)
    f_vec,emin_vec,emax_vec=theAnalysisManager.GetSecondaryElectronSpectrum(root_file,100,0.001*unit.keV, 100.*unit.MeV,
                                logE=True, tree_file_name="secondary_tuple_interaction",
                                weight_correction=weight_correction,extra_cut=extra_cut,
                                prim_flux_emin=prim_flux_emin,
                                prim_flux_emax=prim_flux_emax) 
    evec=0.5*(emin_vec+emax_vec)
    print f_vec
    pl.loglog(evec,f_vec*np.log(eprim_max/eprim_min)/10./nprim,symbol,label=label_name)
    os.chdir(loc_dir)
    
def analysis_de_photo_elec_spectra(test_name,label_name,eprim_min,eprim_max,f_expression="1./ekin_prim",nprim=1.e6,extra_cut="",
                                 prim_flux_emin=0.001,prim_flux_emax=10000000.,symbol=""):
    loc_dir=os.getcwd()
    file_dir=os.path.dirname(__file__)
    if (file_dir != ""):
        os.chdir(file_dir)
    os.chdir("../test_results")
    root_file="%s/%s.root" %(test_name,test_name)
    weight_correction="(ekin_prim)*(%s)" %( f_expression)
    f_vec,emin_vec,emax_vec=theAnalysisManager.GetDeSecondaryElectronSpectrum(root_file,100,0.001*unit.keV, 100.*unit.MeV,
                                logE=True, tree_file_name="secondary_tuple_interaction",
                                weight_correction=weight_correction,extra_cut=extra_cut,
                                prim_flux_emin=prim_flux_emin,
                                prim_flux_emax=prim_flux_emax) 
    evec=0.5*(emin_vec+emax_vec)
    print f_vec
    pl.loglog(evec,f_vec*np.log(eprim_max/eprim_min)/10./nprim,symbol,label=label_name)
    os.chdir(loc_dir)
    
    
    
def analysis_gamma_spectra(test_name,label_name,eprim_min,eprim_max,f_expression="1./ekin_prim",nprim=1.e6,extra_cut="",
                                 prim_flux_emin=0.001,prim_flux_emax=10000000.,symbol=""):
    loc_dir=os.getcwd()
    file_dir=os.path.dirname(__file__)
    if (file_dir != ""):
        os.chdir(file_dir)
    os.chdir("../test_results")
    root_file="%s/%s.root" %(test_name,test_name)
    weight_correction="(ekin_prim)*(%s)" %( f_expression)
    f_vec,emin_vec,emax_vec=theAnalysisManager.GetSecondaryGammaSpectrum(root_file,100,0.001*unit.keV, 100.*unit.MeV,
                                logE=True, tree_file_name="secondary_tuple_interaction",
                                weight_correction=weight_correction,
                                prim_flux_emin=prim_flux_emin,
                                prim_flux_emax=prim_flux_emax) 
    evec=0.5*(emin_vec+emax_vec)
    print f_vec
    pl.loglog(evec,f_vec*np.log(eprim_max/eprim_min)/10./nprim,symbol,label=label_name)
    os.chdir(loc_dir)

    

def analysis_brem_adj_spectra_old(test_names,labels,eprim=10.,nprim=1e5):
    loc_dir=os.getcwd()
    file_dir=os.path.dirname(__file__)
    if (file_dir != ""):
        os.chdir(file_dir)
    os.chdir("../test_results")
    i=0
    for test_name in test_names:
        root_file="%s/%s.root" %(test_name,test_name)
        f_vec,emin_vec,emax_vec=theAnalysisManager.GetAdjSecondaryElecSpectrumForBrem(root_file,100,0.001*unit.keV,
                                                     100.*unit.MeV,
                                logE=True,tree_file_name="secondary_tuple_interaction",eselect_min=eprim*0.95,
                                eselect_max=eprim)
        dE=eprim*0.1
        print f_vec,emin_vec,emax_vec
        evec=0.5*(emin_vec+emax_vec)
        
        pl.loglog(evec,f_vec*np.log(1000)/10/nprim/0.7,label=labels[i])
        i+=1
    #pl.legend(loc="lower left")
    #pl.show()
    os.chdir(loc_dir)
    

def analysis_phot_adj_spectra(test_name,label,eprim_min,eprim_max,f_expression="1./ekin",nprim=1.e6,prim_flux_emin=0.001,
                              prim_flux_emax=10.,symbol="-"):
    loc_dir=os.getcwd()
    file_dir=os.path.dirname(__file__)
    if (file_dir != ""):
        os.chdir(file_dir)
    os.chdir("../test_results")
    i=0
    root_file="%s/%s.root" %(test_name,test_name)
    factor=np.log(eprim_max/eprim_min)
    weight_correction="ekin_prim*%s" %(f_expression)
    f_vec,emin_vec,emax_vec=theAnalysisManager.GetAdjSecondaryElecSpectrumForPhot(root_file,100,0.001*unit.keV,
                                                     100.*unit.MeV,
                                logE=True,tree_file_name="secondary_tuple_interaction",weight_correction=weight_correction,
                                prim_flux_emin=prim_flux_emin,prim_flux_emax=prim_flux_emax)
       
    evec=0.5*(emin_vec+emax_vec)
    print f_vec
    pl.loglog(evec,f_vec*np.log(eprim_max/eprim_min)/10/nprim,symbol,label=label)
    #pl.legend(loc="lower left")
    #pl.show()
    os.chdir(loc_dir)
def analysis_phot_adj_weight_spectra(test_name,label,eprim_min,eprim_max,f_expression="1./ekin",nprim=1.e6,prim_flux_emin=0.001,
                              prim_flux_emax=10.,symbol="-"):
    loc_dir=os.getcwd()
    file_dir=os.path.dirname(__file__)
    if (file_dir != ""):
        os.chdir(file_dir)
    os.chdir("../test_results")
    i=0
    root_file="%s/%s.root" %(test_name,test_name)
    factor=np.log(eprim_max/eprim_min)
    weight_correction="ekin_prim*%s" %(f_expression)
    f_vec,wmin_vec,wmax_vec=theAnalysisManager.GetAdjSecondaryWeightSpectrumForPhot(root_file,100,0.01,
                                                     1000.,
                                logW=True,tree_file_name="secondary_tuple_interaction",weight_correction=weight_correction,
                                prim_flux_emin=prim_flux_emin,prim_flux_emax=prim_flux_emax)
       
    wvec=0.5*(wmin_vec+wmax_vec)
    print f_vec
    pl.loglog(wvec,f_vec*np.log(eprim_max/eprim_min)/10/nprim,symbol,label=label)
    #pl.legend(loc="lower left")
    #pl.show()
    os.chdir(loc_dir)




def analysis_de_phot_adj_spectra(test_name,label,eprim_min,eprim_max,f_expression="1./ekin",nprim=1.e6,prim_flux_emin=0.001,
                              prim_flux_emax=10.,symbol="-"):
    loc_dir=os.getcwd()
    file_dir=os.path.dirname(__file__)
    if (file_dir != ""):
        os.chdir(file_dir)
    os.chdir("../test_results")
    i=0
    root_file="%s/%s.root" %(test_name,test_name)
    factor=np.log(eprim_max/eprim_min)
    weight_correction="ekin_prim*%s" %(f_expression)
    f_vec,emin_vec,emax_vec=theAnalysisManager.GetAdjSecondaryDeSpectrumForPhot(root_file,100,0.001*unit.keV,
                                                     100.*unit.MeV,
                                logE=True,tree_file_name="secondary_tuple_interaction",weight_correction=weight_correction,
                                prim_flux_emin=prim_flux_emin,prim_flux_emax=prim_flux_emax)
       
    evec=0.5*(emin_vec+emax_vec)
    pl.loglog(evec,f_vec*np.log(eprim_max/eprim_min)/10/nprim,symbol,label=label)
    #pl.legend(loc="lower left")
    #pl.show()
    os.chdir(loc_dir)

def analysis_brem_adj_spectra(test_name,label,eprim_adj_min,eprim_adj_max,f_expression="1./ekin",prim_flux_emin=0.001,
                              prim_flux_emax=10.,nprim=1.e6,extra_cut="",symbol=""):
    loc_dir=os.getcwd()
    file_dir=os.path.dirname(__file__)
    if (file_dir != ""):
        os.chdir(file_dir)
    os.chdir("../test_results")
    i=0
    root_file="%s/%s.root" %(test_name,test_name)
    factor=np.log(eprim_adj_max/eprim_adj_min)
    weight_correction="ekin_prim*%s" %(f_expression)
    f_vec,emin_vec,emax_vec=theAnalysisManager.GetAdjSecondaryElecSpectrumForBrem(root_file,100,0.001*unit.keV,
                                                     100.*unit.MeV,
                                logE=True,tree_file_name="secondary_tuple_interaction",
                                weight_correction=weight_correction,extra_cut=extra_cut,
                                prim_flux_emin=prim_flux_emin,prim_flux_emax=prim_flux_emax)
       
    evec=0.5*(emin_vec+emax_vec)
    print f_vec
    pl.loglog(evec,f_vec*np.log(eprim_adj_max/eprim_adj_min)/10/nprim,symbol,label=label)
    #pl.legend(loc="lower left")
    #pl.show()
    os.chdir(loc_dir)


def analyse_brem():
    emin=.0001*unit.MeV
    emax=1000.*unit.MeV
    pl.figure()
    flux_exp1="1."
    flux_exp2="1."
    prim_flux_emin=1.
    prim_flux_emax=10.
    analysis_elec_spectra("Ebrem_0.0001_1000MeV_Tantalum","Ebrem",emin,emax,f_expression=flux_exp1,nprim=1.e5,prim_flux_emin=prim_flux_emin,prim_flux_emax=prim_flux_emax,symbol="k")
    analysis_brem_adj_spectra("InvEbrem_0.0001_1000MeV_Tantalum","InvEbrem",emin,emax,f_expression=flux_exp2,nprim=1.e5,prim_flux_emin=prim_flux_emin,prim_flux_emax=prim_flux_emax,symbol="k^")
    pl.xlabel("kinetic energy")
    pl.ylabel("secondary gamma Flux [a.u]")
    pl.legend(loc="lower left")
    pl.savefig("brem_test.png")
    pl.savefig("brem_test.pdf")
"""
analyse_brem()
pl.show()
    
    
"""

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








#analysis_phot_adj_spectra("InvPEEffect","InvPEEffect",.01*unit.MeV,10.*unit.MeV,f_expression="1.",nprim=1.e6)
flux_exp11="exp(-5.*ekin_prim)"
flux_exp22="exp(-5.*ekin)"
flux_exp1="1."
flux_exp2="1."

"""
pl.figure()
prim_flux_emin=.01
prim_flux_emax=.1
analysis_ephoto_elec_spectra("Compton_0.0001_1000MeV_Tantalum","Fwd 10-100 keV",.001,10.,f_expression=flux_exp1,nprim=1.e6,prim_flux_emin=prim_flux_emin,prim_flux_emax=prim_flux_emax,symbol="k")
analysis_phot_adj_spectra("InvCompt1_0.0001_1000MeV_Tantalum","Adjoint 10-100 keV",.001,10.,f_expression=flux_exp2,nprim=1.e6,prim_flux_emin=prim_flux_emin,prim_flux_emax=prim_flux_emax,symbol="k^")

prim_flux_emin=.1
prim_flux_emax=1.
analysis_ephoto_elec_spectra("Compton_0.0001_1000MeV_Tantalum","Fwd 100keV-1 MeV",.001,10.,f_expression=flux_exp1,nprim=1.e6,prim_flux_emin=prim_flux_emin,prim_flux_emax=prim_flux_emax,symbol="b")
analysis_phot_adj_spectra("InvCompt1_0.0001_1000MeV_Tantalum","Adjoint 100keV-1 MeV",.001,10.,f_expression=flux_exp2,nprim=1.e6,prim_flux_emin=prim_flux_emin,prim_flux_emax=prim_flux_emax,symbol="b^")

prim_flux_emin=1.
prim_flux_emax=10.
analysis_ephoto_elec_spectra("Compton_0.0001_1000MeV_Tantalum","Fwd 1-10 MeV",.001,10.,f_expression=flux_exp1,nprim=1.e6,prim_flux_emin=prim_flux_emin,prim_flux_emax=prim_flux_emax,symbol="r")
analysis_phot_adj_spectra("InvCompt1_0.0001_1000MeV_Tantalum","Adjoint 1-10 MeV",.001,10.,f_expression=flux_exp2,nprim=1.e6,prim_flux_emin=prim_flux_emin,prim_flux_emax=prim_flux_emax,symbol="r^")

prim_flux_emin=10.
prim_flux_emax=100.
analysis_ephoto_elec_spectra("Compton_0.0001_1000MeV_Tantalum","Forward 10-100 MeV",.001,10.,f_expression=flux_exp1,nprim=1.e6,prim_flux_emin=prim_flux_emin,prim_flux_emax=prim_flux_emax,symbol="k")
analysis_phot_adj_spectra("InvCompt1_0.0001_1000MeV_Tantalum","Adjoint 10-100 MeV",.001,10.,f_expression=flux_exp2,nprim=1.e6,prim_flux_emin=prim_flux_emin,prim_flux_emax=prim_flux_emax,symbol="k^")
pl.xlabel("kinetic energy")
pl.ylabel("secondary electron Flux [a.u]")
pl.legend(loc="lower left")
pl.savefig("compton_test.png")
pl.savefig("compton_test.pdf")



pl.figure()

prim_flux_emin=.01
prim_flux_emax=.1
analysis_elec_spectra("PEEffect_0.0001_1000MeV_Tantalum","Fwd 10-100 keV",.0001,1000.,f_expression=flux_exp1,nprim=1.e6,prim_flux_emin=prim_flux_emin,prim_flux_emax=prim_flux_emax,symbol="k")
analysis_phot_adj_spectra("InvPEEffect_0.0001_1000MeV_Tantalum","Adjoint 10-100 keV",.0001,1000.,f_expression=flux_exp2,nprim=1.e6,prim_flux_emin=prim_flux_emin,prim_flux_emax=prim_flux_emax,symbol="k^")


prim_flux_emin=.1
prim_flux_emax=1.
analysis_elec_spectra("PEEffect_0.0001_1000MeV_Tantalum","Fwd 100keV-1 MeV",.0001,1000.,f_expression=flux_exp1,nprim=1.e6,prim_flux_emin=prim_flux_emin,prim_flux_emax=prim_flux_emax,symbol="b")
analysis_phot_adj_spectra("InvPEEffect_0.0001_1000MeV_Tantalum","Adjoint 100keV-1 MeV",.0001,1000.,f_expression=flux_exp2,nprim=1.e6,prim_flux_emin=prim_flux_emin,prim_flux_emax=prim_flux_emax,symbol="b^")

prim_flux_emin=1.
prim_flux_emax=10.
analysis_elec_spectra("PEEffect_0.0001_1000MeV_Tantalum","Fwd 1-10 MeV",.0001,1000.,f_expression=flux_exp1,nprim=1.e6,prim_flux_emin=prim_flux_emin,prim_flux_emax=prim_flux_emax,symbol="r")
analysis_phot_adj_spectra("InvPEEffect_0.0001_1000MeV_Tantalum","Adjoint 1-10 MeV",.0001,1000.,f_expression=flux_exp2,nprim=1.e6,prim_flux_emin=prim_flux_emin,prim_flux_emax=prim_flux_emax,symbol="r^")

prim_flux_emin=10.
prim_flux_emax=100.
analysis_elec_spectra("PEEffect_0.0001_1000MeV_Tantalum","Fwd 10-100 MeV",.0001,1000.,f_expression=flux_exp1,nprim=1.e6,prim_flux_emin=prim_flux_emin,prim_flux_emax=prim_flux_emax,symbol="k")
analysis_phot_adj_spectra("InvPEEffect_0.0001_1000MeV_Tantalum","Adjoint 10-100 MeV",.0001,1000.,f_expression=flux_exp2,nprim=1.e6,prim_flux_emin=prim_flux_emin,prim_flux_emax=prim_flux_emax,symbol="k^")

pl.xlabel("kinetic energy")
pl.ylabel("secondary electron Flux [a.u]")
pl.legend(loc="lower left")
pl.savefig("peffect_test.png")
pl.savefig("peffect_test.pdf")


pl.figure()

prim_flux_emin=.01
prim_flux_emax=.1
analysis_de_photo_elec_spectra("PEEffect_0.0001_1000MeV_Tantalum","Fwd 10-100 keV",.0001,1000.,f_expression=flux_exp1,nprim=1.e6,prim_flux_emin=prim_flux_emin,prim_flux_emax=prim_flux_emax,symbol="k")
analysis_phot_adj_spectra("InvPEEffect_0.0001_1000MeV_Tantalum","Adjoint 10-100 keV",.0001,1000.,f_expression=flux_exp2,nprim=1.e6,prim_flux_emin=prim_flux_emin,prim_flux_emax=prim_flux_emax,symbol="k^")


prim_flux_emin=.1
prim_flux_emax=1.
analysis_de_photo_elec_spectra("PEEffect_0.0001_1000MeV_Tantalum","Fwd 100keV-1 MeV",.0001,1000.,f_expression=flux_exp1,nprim=1.e6,prim_flux_emin=prim_flux_emin,prim_flux_emax=prim_flux_emax,symbol="b")
analysis_de_phot_adj_spectra("InvPEEffect_0.0001_1000MeV_Tantalum","Adjoint 100keV-1 MeV",.0001,1000.,f_expression=flux_exp2,nprim=1.e6,prim_flux_emin=prim_flux_emin,prim_flux_emax=prim_flux_emax,symbol="b^")

prim_flux_emin=1.
prim_flux_emax=10.
analysis_de_photo_elec_spectra("PEEffect_0.0001_1000MeV_Tantalum","Fwd 1-10 MeV",.0001,1000.,f_expression=flux_exp1,nprim=1.e6,prim_flux_emin=prim_flux_emin,prim_flux_emax=prim_flux_emax,symbol="r")
analysis_de_phot_adj_spectra("InvPEEffect_0.0001_1000MeV_Tantalum","Adjoint 1-10 MeV",.0001,1000.,f_expression=flux_exp2,nprim=1.e6,prim_flux_emin=prim_flux_emin,prim_flux_emax=prim_flux_emax,symbol="r^")

prim_flux_emin=10.
prim_flux_emax=100.
analysis_de_photo_elec_spectra("PEEffect_0.0001_1000MeV_Tantalum","Fwd 10-100 MeV",.0001,1000.,f_expression=flux_exp1,nprim=1.e6,prim_flux_emin=prim_flux_emin,prim_flux_emax=prim_flux_emax,symbol="k")
analysis_de_phot_adj_spectra("InvPEEffect_0.0001_1000MeV_Tantalum","Adjoint 10-100 MeV",.0001,1000.,f_expression=flux_exp2,nprim=1.e6,prim_flux_emin=prim_flux_emin,prim_flux_emax=prim_flux_emax,symbol="k^")

pl.xlabel("kinetic energy")
pl.ylabel("secondary electron Flux [a.u]")
pl.legend(loc="lower left")
pl.savefig("de_peffect_test.png")
pl.savefig("de_peffect_test.pdf")

pl.figure()
"""
prim_flux_emin=.38
prim_flux_emax=.39
analysis_phot_adj_weight_spectra("InvPEEffect_0.0001_1000MeV_Tantalum","Adjoint ",.0001,1000.,f_expression=flux_exp2,nprim=1.e6,prim_flux_emin=prim_flux_emin,prim_flux_emax=prim_flux_emax,symbol="k^")


prim_flux_emin=1.6
prim_flux_emax=1.7
analysis_phot_adj_weight_spectra("InvPEEffect_0.0001_1000MeV_Tantalum","Adjoint ",.0001,1000.,f_expression=flux_exp2,nprim=1.e6,prim_flux_emin=prim_flux_emin,prim_flux_emax=prim_flux_emax,symbol="r^")


pl.show()

pl.xlabel("kinetic energy")
pl.ylabel("weight [a.u]")
pl.legend(loc="lower left")
pl.savefig("adjoint_weight_peffect_test.png")
pl.savefig("adjoint_weight_peffect_test.pdf")





"""
pl.figure()
prim_flux_emin=.01
prim_flux_emax=1.

analysis_ephoto_elec_spectra("Eioni_0.0001_1000MeV_Tantalum","Fwd prim e-",.0001,1000.,f_expression=flux_exp1,nprim=1.e6,extra_cut=" && ekin > ekin_prim/2.",
                                                                    prim_flux_emin=prim_flux_emin,prim_flux_emax=prim_flux_emax,symbol="k")
analysis_ephoto_elec_spectra("Eioni_0.0001_1000MeV_Tantalum","Fwd second e-",.0001,1000.,f_expression=flux_exp1,nprim=1.e6,extra_cut=" && ekin < ekin_prim/2.",
                                                                    prim_flux_emin=prim_flux_emin,prim_flux_emax=prim_flux_emax,symbol="r")


#analysis_ephoto_elec_spectra("Eioni_10_100keV","Eioni_10_100keV",.1,1.,f_expression=flux_exp1,nprim=1.e6,extra_cut=" && ekin < ekin_prim/2.")
#analysis_ephoto_elec_spectra("Eioni","Eioni",emin,emax,f_expression=flux_exp1,nprim=1.e5,extra_cut=" && ekin < ekin_prim/2.")


#analysis_brem_adj_spectra("InvEioni_10_100keV","InvEioni_10_100keV",.01,.1,f_expression=flux_exp1,nprim=1.e6,prim_flux_emin=prim_flux_emin,prim_flux_emax=prim_flux_emax)
#analysis_brem_adj_spectra("InvEioni1_10_100keV","InvEioni10_100keV",.01,.1,f_expression=flux_exp1,nprim=1.e6,prim_flux_emin=prim_flux_emin,prim_flux_emax=prim_flux_emax)

analysis_brem_adj_spectra("InvEioni_0.0001_1000MeV_Tantalum","Adjoint prim e-",0.0001,1000,f_expression=flux_exp2,nprim=1.e6,
                                                                                                prim_flux_emin=prim_flux_emin,
                                                                                                prim_flux_emax=prim_flux_emax,symbol="k^")
analysis_brem_adj_spectra("InvEioni1_0.0001_1000MeV_Tantalum","Adjoint second e-",0.0001,1000,f_expression=flux_exp2,nprim=1.e6,
                                                    prim_flux_emin=prim_flux_emin,prim_flux_emax=prim_flux_emax,symbol="r^")

pl.xlabel("kinetic energy")
pl.ylabel("electron Flux [a.u]")
pl.legend(loc="lower left")
pl.savefig("ionisation_test.png")
pl.savefig("ionisation_test.pdf")

"""


