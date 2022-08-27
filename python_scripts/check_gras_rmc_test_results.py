#!/usr/bin/python
import sys,os
#It is important to set the directory, where  the subdirectory python_utilities is located 
# in the python  path thi is done below
sys.path.append("..")
from python_utilities import GRASTestor
from python_utilities.Plotters import Plotter
from python_utilities import SpenvisCSVFileHandler
import numpy as np
import pylab as pl
from python_utilities import hepunit as unit
import scipy as scp

from ROOT import TF1,gRandom, gROOT,gPad, gStyle, gDirectory,TH2F,TH3F,\
                    TClass,TEntryList, TNtuple,TObject,TDirectory, TObjLink,\
                    TCanvas, TPad, TH1F , TFile, TLegend, TGraph, TVector,\
                    Math, TGraph,Double,Long,TEventList



#Create the SpernvisCSVFileHandler instance
###########################
theSpenvisCSVFileHandler=SpenvisCSVFileHandler.SpenvisCSVFileHandler()
thePlotter=Plotter()


def MakeArrayForBranchAndSetAdress(root_tree,branch_name):
        vector=None
        type_str=root_tree.GetBranch(branch_name).GetTitle()[-2:].upper()
        if (type_str=="/F"):
            vector=np.array( [ 0. ],'f' )
        elif (type_str=="/I"):
            vector=np.array( [ 0. ],'i' )
        if vector is not None:
            root_tree.SetBranchAddress(branch_name,vector)
        return vector

def check_adj_dose_from_tuple(test_name,f_diff_func=None,sensitive_volume="GreatBoxBoard1Component",
                              nb_event=10e4):
    script_dir=os.path.dirname(__file__)
    cwd=os.getcwd()
    file_dir=os.path.dirname(__file__)
    if (file_dir != ""):
        os.chdir(file_dir) 
    os.chdir("..")
    test_dir="%s/test_results/%s" %(os.getcwd(),test_name)
    RootFileName="%s/%s_adj_run.root" %(test_dir,test_name)
    root_file=TFile(RootFileName,"READ")
    adjoint_source_tuple=root_file.Get("tuple_adjoint_prim_on_external_source")
    map_adjoint_tuple=root_file.Get("adjoint_EdepIn%s_tuple_dose" %(sensitive_volume))
    dose_tuple=root_file.Get("EdepIn%s_tuple_dose" %(sensitive_volume))
    
    
    n_dose_tuple=dose_tuple.GetEntries()
    dose_v= MakeArrayForBranchAndSetAdress(dose_tuple,"dose")
    primarykine_v= MakeArrayForBranchAndSetAdress(dose_tuple,"primarykine")
    weight_v= MakeArrayForBranchAndSetAdress(dose_tuple,"weight")
    #pdg_v= MakeArrayForBranchAndSetAdress(dose_tuple,"pdg")
    adj_id_v= MakeArrayForBranchAndSetAdress(map_adjoint_tuple,"adj_prim_nb")
    
    i_adj_source_tuple=0
    n_adj_source_tuple=adjoint_source_tuple.GetEntries()
    pdg0_v=MakeArrayForBranchAndSetAdress(adjoint_source_tuple,"pdg0")
    pdg1_v=MakeArrayForBranchAndSetAdress(adjoint_source_tuple,"pdg1")
    kine0_v=MakeArrayForBranchAndSetAdress(adjoint_source_tuple,"kine0")
    kine1_v=MakeArrayForBranchAndSetAdress(adjoint_source_tuple,"kine1")
    weight1_v=MakeArrayForBranchAndSetAdress(adjoint_source_tuple,"weight")
    evt_v=MakeArrayForBranchAndSetAdress(adjoint_source_tuple,"evt")
    dose_vec=[]
    d_dose_vec=[]
    evt_vec=[]
    icount=0
    dose=0.
    n_adj_source_tuple=10000
    for i in range(n_adj_source_tuple):
        dose_tuple.GetEntry(i)
        map_adjoint_tuple.GetEntry(i)
        adjoint_source_tuple.GetEntry(adj_id_v[0]-1)
        if pdg1_v[0] ==11:
            icount+=1
            """
            print "Dose , energy",dose_v[0],primarykine_v[0],weight_v[0]
            print pdg1_v[0],kine1_v[0], pdg0_v[0],kine0_v[0],weight1_v[0]
            """
            d_dose=dose_v[0]*weight1_v[0]*f_diff_func(kine1_v[0])
            dose+=d_dose
            d_dose_vec+=[d_dose]
            dose_vec+=[dose/evt_v[0]]
            evt_vec+=[evt_v[0]]
    
    pl.subplot(3,1,1)
    pl.plot(np.arange(len(dose_vec)),dose_vec)
    indices=np.argsort(d_dose_vec)
    sort_d_dose_vec=np.sort(d_dose_vec)
    dose_vec1=[]
    dose1=0.
    for i in range(len(d_dose_vec)):
        dose1+=sort_d_dose_vec[i]
        dose_vec1+=[dose1/evt_vec[i]]
    
    #pl.plot(np.array(evt_vec)+1.,dose_vec1)
    pl.plot(np.arange(len(dose_vec)-3),dose_vec1[0:-3])
        
        
    
   
    
    #pl.ylim(np.max(dose_vec)*np.array([0.00001,3.]))
    
    
     #Plot the conv curve
    table_name="%s/%s_adj_runConvergence_EdepIn%s.txt" %(test_dir,test_name,sensitive_volume)
    
    res=np.loadtxt(table_name,skiprows=1)
    
    pl.subplot(3,1,2)
    pl.loglog(res[:,-2],res[:,-3])
    
    print res[:,-2],res[:,-3]
    pl.ylabel("precision [%]" )
    pl.xlabel("computing time")
    pl.subplot(3,1,3)
    pl.semilogx(res[:,-2],res[:,0])
    pl.semilogx(res[:,-2],(100.+res[:,-3])*res[-1,0]/100.)
    pl.semilogx(res[:,-2],(100.-res[:,-3])*res[-1,0]/100.)
    pl.ylim([0.5*res[-1,0],1.5*res[-1,0]])
    
    pl.ylabel("dose [a.u.]" )
    pl.xlabel("computing time")
    
    
    
    pl.show()
    os.chdir(cwd)
            
def fdiff_exp(ekin):
    return np.exp(-ekin/1.)

def fdiff_juice(ekin):
    
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
             
    return np.exp(scp.interp(ekin,evec1,np.log(fdiff_vec1)))
    



#check_adj_dose_from_tuple("test_ngrm_juice_2mm_Ta_housing",f_diff_func=fdiff_juice,sensitive_volume="GreatBoxBoard1Component")
#check_adj_dose_from_tuple("test_ngrm_juice_2mm_Ta_housing_no_adj_PE_and_Compton",f_diff_func=fdiff_juice,sensitive_volume="GreatBoxBoard1Component")
#check_adj_dose_from_tuple("test_ngrm_juice_2mm_Ta_housing_fwdCS_mode_new",f_diff_func=fdiff_juice,sensitive_volume="GreatBoxBoard1Component")
            
#check_adj_dose_from_tuple("test_ngrm_juice_2mm_Ta_housing_fwdCS_mode_long",f_diff_func=fdiff_juice,sensitive_volume="GreatBoxBoard1Component")
            
             
            
        
        
        
        
    
    
    
    

def check_results(test_name,show_plots=True,save_figure=True,
                  sensitive_layer_index=1,sensitive_volume=None,
                  suffix_adj_run="",suffix_fwd_run=""):
    """
    Will compare all histos from adjoint and forward sim 
    """
    script_dir=os.path.dirname(__file__)
    cwd=os.getcwd()
    file_dir=os.path.dirname(__file__)
    if (file_dir != ""):
        os.chdir(file_dir) 
    os.chdir("..")
    test_dir="%s/test_results/%s" %(os.getcwd(),test_name)
    os.chdir(cwd)
    adjoint_csv_file="%s/%s_adj_run%s_Spectrum1.csv" %(test_dir,test_name,suffix_adj_run)
    print adjoint_csv_file
    fwd_csv_file="%s/%s_fwd_run%s.csv" %(test_dir,test_name,suffix_fwd_run)
    
    
    if not (os.path.exists(adjoint_csv_file) and os.path.exists(fwd_csv_file)) :
        print "Test %s does not exists" %(test_name)
        print os.path.exists(adjoint_csv_file),adjoint_csv_file
        print os.path.exists(fwd_csv_file),fwd_csv_file
        return
    
    #Get the test results
    
    ModuleDescription, StatDoubleTable_adj, Histo1Ds_adj=theSpenvisCSVFileHandler.ReadGRASFile(adjoint_csv_file)
    ModuleDescription, StatDoubleTable_fwd, Histo1Ds_fwd=theSpenvisCSVFileHandler.ReadGRASFile(fwd_csv_file)
    adj_dose=0.
    fwd_dose=0.
    if sensitive_volume is None:
        adj_dose=StatDoubleTable_adj["EdepInLayer%i" %(sensitive_layer_index+1)]["value"]
        fwd_dose=StatDoubleTable_fwd["EdepInLayer%i" %(sensitive_layer_index+1)]["value"]
    else:
        adj_dose=StatDoubleTable_adj["EdepIn%s" %(sensitive_volume)]["value"]
        fwd_dose=StatDoubleTable_fwd["EdepIn%s" %(sensitive_volume)]["value"]
    
    
    
    
    #Plot the histos
    
    
    for key in Histo1Ds_fwd.keys():
        for title  in Histo1Ds_fwd[key]:
            fig_base_name=test_dir+"/"+key+"_"+title
            histos=[Histo1Ds_adj[key][title],Histo1Ds_fwd[key][title]]
            thePlotter.PlotListGRAS1DHisto(histos,fig_base_name,1,1,
                                                show_plots=show_plots,titles=title,
                                                all_on_same_plot=True,ylogs=True,
                                                xlogs=True,labels=["adj %.2e MeV" %(adj_dose),"fwd %.2e MeV" %(fwd_dose)])
    
    #Plot the conv curve
    table_name="%s/%s_adj_run%sConvergence_EdepInLayer%i.txt" %(test_dir,test_name,suffix_adj_run,sensitive_layer_index+1)
    if sensitive_volume is not None:
        table_name="%s/%s_adj_run%sConvergence_EdepIn%s.txt" %(test_dir,test_name,suffix_adj_run,sensitive_volume)
    
    res=np.loadtxt(table_name,skiprows=1)
    if save_figure:
        pl.figure(figsize=(8,6))
    pl.subplot(2,1,1)
    pl.loglog(res[:,-2],res[:,-3])
    print res[:,-2],res[:,-3]
    pl.ylabel("precision [%]" )
    pl.xlabel("computing time")
    pl.subplot(2,1,2)
    pl.semilogx(res[:,-2],res[:,0])
    pl.semilogx(res[:,-2],(100.+res[:,-3])*res[-1,0]/100.)
    pl.semilogx(res[:,-2],(100.-res[:,-3])*res[-1,0]/100.)
    pl.ylim([0.5*res[-1,0],1.5*res[-1,0]])
    
    pl.ylabel("dose [a.u.]" )
    pl.xlabel("computing time")
    
    
    """
    pl.subplot(2,1,2)
    pl.plot(res[:,-2],res[:,0]/1000)
    pl.ylabel("edep [MeV]" )
    pl.xlabel("computing time")
    """
    fig_file="%s/%s_adj_runConvergence_EdepInLayer%i.png" %(test_dir,test_name,sensitive_layer_index+1)
    if save_figure:
        pl.savefig(fig_file)
    pl.show()
    

def compare_group_convergence_curves(test_name,show_plots=True,save_figure=True,
                  sensitive_layer_index=1,sensitive_volume=None,
                  suffix_adj_run_tags=[],suffix_fwd_run="",labels=[],styles=[],
                  xmax=1.e4,plot_name="test_convergence"):
    
    
    """
    Will compare all histos from adjoint and forward sim 
    """
    script_dir=os.path.dirname(__file__)
    cwd=os.getcwd()
    file_dir=os.path.dirname(__file__)
    if (file_dir != ""):
        os.chdir(file_dir) 
    os.chdir("..")
    test_dir="%s/test_results/%s" %(os.getcwd(),test_name)
    os.chdir(cwd)
    
    fwd_csv_file="%s/%s_fwd_run%s.csv" %(test_dir,test_name,suffix_fwd_run)
    
    ModuleDescription, StatDoubleTable, Histo1Ds=theSpenvisCSVFileHandler.ReadGRASFile(fwd_csv_file)
    
    fwd_dose=0.
    fwd_dose_error=0.
    if sensitive_volume is None:
        fwd_dose=StatDoubleTable["EdepInLayer%i" %(sensitive_layer_index+1)]["value"]
        fwd_dose_error=StatDoubleTable["EdepInLayer%i" %(sensitive_layer_index+1)]["error"]
    else:
        fwd_dose=StatDoubleTable["EdepIn%s" %(sensitive_volume)]["value"]
        fwd_dose_error=StatDoubleTable["EdepIn%s" %(sensitive_volume)]["error"]
    
    
    file_list=os.listdir(test_dir)
    table_name_mat=[]
    for tag in suffix_adj_run_tags:
        table_name_mat+=[[]]
        for fname in file_list :
            if tag in fname and "Convergence" in fname:
                table_name_mat[-1]+=[fname]

    for i in range(len(table_name_mat)):
        table_name=table_name_mat[i][0]
        res=np.loadtxt("%s/%s" %(test_dir,table_name),skiprows=1)
    
        pl.subplot(2,1,1)
        pl.loglog(res[:,-2],res[:,-3],styles[i],label=labels[i])
        pl.ylabel("precision [%]" )
        pl.xlabel("computing time")
        pl.subplot(2,1,2)
        pl.semilogx(res[:,-2],res[:,0],styles[i])
        pl.ylim([0.5*res[-1,0],1.5*res[-1,0]])
    
        pl.ylabel("dose [a.u.]" )
        pl.xlabel("computing time")
        for table_name in table_name_mat[i][1:]:
            res=np.loadtxt("%s/%s" %(test_dir,table_name),skiprows=1)
    
            pl.subplot(2,1,1)
            pl.loglog(res[:,-2],res[:,-3],styles[i])
            pl.subplot(2,1,2)
            pl.semilogx(res[:,-2],res[:,0],styles[i])
        pl.subplot(2,1,1)
        pl.legend()
        pl.xlim([0.5,xmax])
        pl.subplot(2,1,2)
        x=[0.5,xmax]
        y=[fwd_dose,fwd_dose]
        pl.plot(x,y,"g")
        pl.plot(x,y-fwd_dose_error,"g--")
        pl.plot(x,y+fwd_dose_error,"g--")
        
        
        pl.xlim([0.5,xmax])
        
        
        
    
    if save_figure:
        pl.savefig("%s.pdf" %(plot_name))
        pl.savefig("%s.png" %(plot_name))
        
    pl.show()
    
    
def compare_dose_distribution(test_name,show_plots=True,save_figure=True,
                  sensitive_layer_index=1,sensitive_volume=None,
                  suffix_adj_run_tags=[],suffix_fwd_run="",labels=[],styles=[]):
            
          
                
    """
    Will compare all histos from adjoint and forward sim 
    """
    script_dir=os.path.dirname(__file__)
    cwd=os.getcwd()
    file_dir=os.path.dirname(__file__)
    if (file_dir != ""):
        os.chdir(file_dir) 
    os.chdir("..")
    test_dir="%s/test_results/%s" %(os.getcwd(),test_name)
    os.chdir(cwd)
    
    file_list=os.listdir(test_dir)
    dose_mat=[]
    for tag in suffix_adj_run_tags:
        dose_mat+=[[]]
        names_vec=[]
        for fname in file_list :
            if tag in fname and "Spectrum1" in fname:
                print fname
                ModuleDescription, StatDoubleTable, Histo1Ds=theSpenvisCSVFileHandler.ReadGRASFile("%s/%s" %(test_dir,fname))
                names_vec+=[fname]
                dose=0.
                if sensitive_volume is None:
                    dose=StatDoubleTable["EdepInLayer%i" %(sensitive_layer_index+1)]["value"]
                else:
                    dose=StatDoubleTable["EdepIn%s" %(sensitive_volume)]["value"]
                dose_mat[-1]+=[dose]
        dose_mat[-1]=np.array(dose_mat[-1])    
        ind=np.argmax(dose_mat[-1])
        
        print np.mean(dose_mat[-1]),np.median(dose_mat[-1]),np.std(dose_mat[-1]),dose_mat[-1][ind],names_vec[ind]
         
    
    print dose_mat
                    
                
    """
    for i in range(len(table_name_mat)):
        table_name=table_name_mat[i][0]
        res=np.loadtxt("%s/%s" %(test_dir,table_name),skiprows=1)          
    """
    """
    #Plot the conv curve
    table_name="%s/%s_adj_run%sConvergence_EdepInLayer%i.txt" %(test_dir,test_name,suffix_adj_run,sensitive_layer_index+1)
    if sensitive_volume is not None:
        table_name="%s/%s_adj_run%sConvergence_EdepIn%s.txt" %(test_dir,test_name,suffix_adj_run,sensitive_volume)
    """





def comp_dose_results():
    test_names=["test_Al_1.00mmAlEq_0.01mm_cut",
                "test_Al_5.00mmAlEq_0.01mm_cut",
                "test_Al_10.00mmAlEq_0.01mm_cut",
                "test_Al_1.00mmAlEq_10mm_cut",
                "test_Al_5.00mmAlEq_10mm_cut",
                "test_Al_10.00mmAlEq_10mm_cut",
                "test_Ta_1.00mmAlEq_0.01mm_cut",
                "test_Ta_5.00mmAlEq_0.01mm_cut",
                "test_Ta_10.00mmAlEq_0.01mm_cut",
                "test_Ta_1.00mmAlEq_10mm_cut",
                "test_Ta_5.00mmAlEq_10mm_cut",
                "test_Ta_10.00mmAlEq_10mm_cut"
                
                
                
                ]
    
    script_dir=os.path.dirname(__file__)
    cwd=os.getcwd()
    file_dir=os.path.dirname(__file__)
    if (file_dir != ""):
        os.chdir(file_dir) 
    labels_vec=[]

    
    adj_dose=[]
    fwd_dose=[]
    err_adj_dose=[]
    err_fwd_dose=[]
    for test_name in test_names:
        os.chdir("..")
        test_dir="%s/test_results/%s" %(os.getcwd(),test_name)
        os.chdir(cwd)
        adjoint_csv_file="%s/%s_adj_run_Spectrum1.csv" %(test_dir,test_name)
        fwd_csv_file="%s/%s_fwd_run.csv" %(test_dir,test_name)
    
        if not (os.path.exists(adjoint_csv_file) and os.path.exists(adjoint_csv_file)) :
            print "Test %s does not exists" %(test_name)
            return
        
        
        sensitive_layer_index=1
    
        #Get the test results
        ModuleDescription, StatDoubleTable_adj, Histo1Ds_adj=theSpenvisCSVFileHandler.ReadGRASFile(adjoint_csv_file)
        ModuleDescription, StatDoubleTable_fwd, Histo1Ds_fwd=theSpenvisCSVFileHandler.ReadGRASFile(fwd_csv_file)
        adj_dose+=[StatDoubleTable_adj["EdepInLayer%i" %(sensitive_layer_index+1)]["value"]]
        fwd_dose+=[StatDoubleTable_fwd["EdepInLayer%i" %(sensitive_layer_index+1)]["value"]]
        err_adj_dose+=[StatDoubleTable_adj["EdepInLayer%i" %(sensitive_layer_index+1)]["error"]]
        err_fwd_dose+=[StatDoubleTable_fwd["EdepInLayer%i" %(sensitive_layer_index+1)]["error"]]
    pl.figure(figsize=[12,8])
    pl.rcParams['xtick.labelsize']=18
    pl.rcParams['ytick.labelsize']=18
    pl.rcParams['axes.labelsize']=18
    pl.rcParams['legend.fontsize']=16
    pl.rcParams['font.weight']='bold'
    
    pl.semilogy([1.,5.,10.],fwd_dose[0:3],"k")
    pl.errorbar([1.,5.,10.],fwd_dose[0:3],yerr=err_fwd_dose[0:3],fmt="k",label="Al shielding Forward, 0.01mm cut")
    pl.errorbar([1.,5.,10.],adj_dose[0:3],yerr=err_adj_dose[0:3],fmt="ok",label="Al shielding Reverse, 0.01mm cut")
    pl.errorbar([1.,5.,10.],fwd_dose[3:6],yerr=err_fwd_dose[3:6],fmt="k--",label="Al shielding Forward, 10mm cut")
    pl.errorbar([1.,5.,10.],adj_dose[3:6],yerr=err_adj_dose[3:6],fmt="*k",label="Al shielding Reverse, 10mm cut")
    pl.errorbar([1.,5.,10.],fwd_dose[6:9],yerr=err_fwd_dose[6:9],fmt="b",label="Ta shielding Forward, 0.01mm cut")
    pl.errorbar([1.,5.,10.],adj_dose[6:9],yerr=err_adj_dose[6:9],fmt="ob",label="Ta shielding Reverse, 0.01mm cut")
    pl.errorbar([1.,5.,10.],fwd_dose[9:],yerr=err_fwd_dose[9:],fmt="b--",label="Ta shielding Forward, 10mm cut")
    pl.errorbar([1.,5.,10.],adj_dose[9:],yerr=err_adj_dose[9:],fmt="*b",label="Ta shielding Reverse, 10mm cut")
    pl.legend()
    pl.xlabel("Al thickness [mm]")
    pl.ylabel("Energy deposited  [MeV]")
    
    
    
    
    
    
    pl.xlim([0.5,12.])
    pl.savefig("spherical_shielding.png")
    pl.show()
    
    
        
    
    
                


def check_results_for_several_sims(test_names,labels_sim,show_plots=True,save_figure=True,
                  sensitive_layer_index=1,prefix_plot="test_plot",fmts=None,
                  sensitive_volume_name=None,suffix_adj_run="3",suffix_fwd_run=""):
    """
    Will compare all histos from adjoint and forward sim 
    """
    script_dir=os.path.dirname(__file__)
    cwd=os.getcwd()
    file_dir=os.path.dirname(__file__)
    if (file_dir != ""):
        os.chdir(file_dir) 
    labels_vec=[]
    for label in labels_sim:
        labels_vec+=[label+" adj",label+" fwd"]
    histos_fwd={}
    histos_adj={}
    for test_name in test_names:
        os.chdir("..")
        test_dir="%s/test_results/%s" %(os.getcwd(),test_name)
        os.chdir(cwd)
        adjoint_csv_file="%s/%s_adj_run%s_Spectrum1.csv" %(test_dir,test_name,suffix_adj_run)
        fwd_csv_file="%s/%s_fwd_run.csv" %(test_dir,test_name)
    
        if not (os.path.exists(adjoint_csv_file) and os.path.exists(adjoint_csv_file)) :
            print "Test %s does not exists" %(test_name)
            return
    
        #Get the test results
        ModuleDescription, StatDoubleTable_adj, Histo1Ds_adj=theSpenvisCSVFileHandler.ReadGRASFile(adjoint_csv_file)
        ModuleDescription, StatDoubleTable_fwd, Histo1Ds_fwd=theSpenvisCSVFileHandler.ReadGRASFile(fwd_csv_file)
        """"
        adj_dose=StatDoubleTable_adj["EdepInLayer%i" %(sensitive_layer_index+1)]["value"]
        fwd_dose=StatDoubleTable_fwd["EdepInLayer%i" %(sensitive_layer_index+1)]["value"]
        err_adj_dose=StatDoubleTable_adj["EdepInLayer%i" %(sensitive_layer_index+1)]["error"]
        fwd_adj_dose=StatDoubleTable_fwd["EdepInLayer%i" %(sensitive_layer_index+1)]["error"]
        print adj_dose,err_adj_dose
        print fwd_dose,err_fwd_dose
        """
        
        
    
        #Plot the histos
        for key in Histo1Ds_fwd.keys():
            if key not in histos_fwd:
                histos_fwd[key]=[]
                histos_adj[key]=[]
            histos_fwd[key]+=[Histo1Ds_fwd[key]]
            histos_adj[key]+=[Histo1Ds_adj[key]]
    i=0
    pl.figure(figsize=[12,6])
    pl.rcParams['xtick.labelsize']=18
    pl.rcParams['ytick.labelsize']=18
    pl.rcParams['axes.labelsize']=18
    pl.rcParams['legend.fontsize']=16
    pl.rcParams['font.weight']='bold'
    for key in histos_fwd.keys():
        for title in histos_fwd[key][0]:
            histos=[]
            i=0
            for test_name in test_names:
                histos+=[histos_adj[key][i][title],histos_fwd[key][i][title]]
                i+=1
            fig_base_name=prefix_plot+"_"+key+"_"+title
            thePlotter.PlotListGRAS1DHisto(histos,fig_base_name,1,1,
                                                show_plots=show_plots,titles=title,
                                                all_on_same_plot=True,ylogs=True,
                                                xlogs=True,labels=labels_vec,
                                                fmts=fmts)
    
   
    
    
    
    
    
    
#check_results("test_1D_5mmVacuum")
#check_results("test_1D_0.5mmAl")  
"""
al_thicknesses=[0.5,1.,5.,10.]
E0_vec=[0.25,0.5,1.,5.,10.,20.]
materials=["Aluminum","Tantalum"]
for al_thickness in al_thicknesses[-2:]:
    for E0 in E0_vec[-2:]:
        for material in materials:
            test_name="test_1D_%.4emmeqAl_%s_%.4eMeVE0" %(al_thickness,material,E0)
            check_results(test_name,show_plots=False) 
            test_name="test_1D_%.4emmeqAl_%s_%.4eMeVE0_big_cut" %(al_thickness,material,E0)
            check_results(test_name,show_plots=False)

#check_results("test_1D_10mm_Al_E020MeV_cut_big",show_plots=False)
"""
"""
al_thicknesses=np.array([0.5,1.,5.,10.])
al_thicknesses=np.array([5.])
materials=["Aluminum","Tantalum"]
#materials=["Vacuum"]
E0_vec=np.array([1.,5.,10.,20.])
emax_vec=np.array([10.,20.,50.,100.,1000.])
#emax_vec=np.array([1.,10.,100.,1000.])
pl.figure(figsize=(6,4))
for al_thickness in al_thicknesses:
    for material in materials[1:]:
        for l in [3]:
            E0=E0_vec[l]
            emax=emax_vec[l]
            tag="noMS_"
            test_name="test_1D_%.4emmeqAl_%s_%se0_%.2f_emax_%.2f_new_along" %(al_thickness,material,tag,E0,emax)
            
            check_results(test_name,show_plots=False)
            
            for tag in ["MS_io_and_brem_","noMS_io_and_brem_","noMS_","","noCompton_","noPEEFFECT_","_"][0:1]:
        
                test_name="test_1D_%.4emmeqAl_%s_%se0_%.2f_emax_%.2f" %(al_thickness,material,tag,emax)
                check_results(test_name,show_plots=False)
            
            for tag in ["","MS_io_and_brem_","noMS_"]:
                test_name="test_1D_%.4emmeqAl_%s_%sJuice_emax_%.2f" %(al_thickness,material,tag,emax)
                check_results(test_name,show_plots=False,save_figure=False)

pl.savefig("convergence.png")            
"""        
"""
check_results("test_Ta_1.00mm_noAdjointComptonPhoto_juice",
              show_plots=True,
              save_figure=True,
              sensitive_layer_index=1)
"""
"""
check_results_for_several_sims(["test_Ta_2.00mm_io_brem_juice",
                               "test_Ta_2.00mm_io_brem_noMS_juice",
                               "test_Ta_2.00mm_juice"],
                            ["io brem ","io brem no MS","all processes"],
              show_plots=True,
              save_figure=True,
              sensitive_layer_index=1,fmts=["k-","k--","r-","r--","b-","b--"])



check_results_for_several_sims(["test_Al_5.00mm_io_brem_E05.00e+01MeV_Emin1.00e+01MeV_Emax1.00e+02MeV"],
                            ["10 MeV 100 MeV "],
              show_plots=True,
              save_figure=True,
              sensitive_layer_index=1,fmts=["k-","k--","r-","r--","b-","b--"])


check_results_for_several_sims(["test_payload_juice_big_s_volume_Ta3mm"],
                            ["juice"],
              show_plots=True,
              save_figure=True,
              fmts=["k-","k--","r-","r--","b-","b--"])
"""

"""
check_results_for_several_sims(["test_Al_5.00mm_io_brem_juice",
                               "test_Al_5.00mm_juice"],
                            ["io brem ","all processes"],
              show_plots=False,
              save_figure=True,prefix_plot="../plots/test_Al_5.00",
              sensitive_layer_index=1,fmts=["k-","k--","r-","r--","b-","b--"],
              
              )

check_results_for_several_sims(["test_Al_5.00mm_io_brem_juice",
                               "test_Al_5.00mm_io_brem_noMS_juice"],
                            ["io brem ","io brem no MS"],
              show_plots=False,
              save_figure=True,prefix_plot="../plots/test_Al_5.00_io_brem",
              sensitive_layer_index=1,fmts=["k-","k--","r-","r--","b-","b--"],
              
              )





check_results("test_Al_0.10mm_io_brem_E01.00e+00MeV_Emin1.00e-02MeV_Emax1.00e+00MeV",
              show_plots=True,
              save_figure=True,
              sensitive_layer_index=1,
              sensitive_volume="component1")

"""


"""
comp_dose_results()
"
"""
test_name="test_ngrm_juice_2mm_Ta_housing_fwdCS_mode_test"
test_name="test_Ta_5.00mmAlEq_0.01mm_cut"
test_name="test_Ta_10.00mmAlEq_0.01mm_cut"
test_name="test_Ta_10.00mmAlEq_0.01mm_cut_NoCompton_PhotoElec"
test_name="test_Ta_10.00mmAlEq_0.01mm_cut_NoPhoto"
test_name="test_Ta_10.00mmAlEq_0.01mm_cut_AllNoSecondaryElectronsThirdNoAdjointPhotoElec"
test_name="test_ngrm_juice_2mm_Ta_housing_NoElectronSecondaryButDelta"
test_name="test_ngrm_juice_2mm_Ta_housing_NoAdjointComptonandPElec"
#test_name="test_ngrm_juice_2mm_Ta_housing_OnlyIonisationAndBrem"
#test_name="test_ngrm_juice_2mm_Ta_housing_NoElectronSecondaryButDelta"
#test_name="test_ngrm_juice_4mm_Ta_housing_NoAdjComptonAndPhoto"
test_name="test_ngrm_juice_15mm_Al_housing_NoAdjComptonAndPhoto"
test_name="test_ngrm_gallileo_0.5mm_Ta_housing_NoAdjointComptonAndPElec"
test_name="test_ngrm_gallileo_2mm_Al_housing_NoAdjointComptonAndPElec"
#test_name="test_ngrm_gallileo_0.5mm_Ta_housing_NoAdjointComptonandPElec"
#test_name="test_ngrm_gallileo_1mm_Ta_housing_NoElectronToGammaReverse"


check_results(test_name,
              show_plots=True,
              save_figure=True,
              sensitive_layer_index=1
              ,sensitive_volume="GreatBoxBoard1Component",
              suffix_adj_run="1e5evts_splitf10_10",suffix_fwd_run="1")


compare_group_convergence_curves(test_name,show_plots=True,save_figure=True,
                  sensitive_layer_index=1,sensitive_volume="GreatBoxBoard1Component",
                  suffix_adj_run_tags=["run1e5evts_splitf1_","run1e5evts_splitf10_"],
                  suffix_fwd_run="1",labels=["fsplit 1","fsplit 10"],styles=["k","r","g"],
                  plot_name="test_convergence_ngrm_gallileo_Ta0.5mm_housing"
                  )

"""
compare_dose_distribution(test_name,show_plots=True,save_figure=True,
                  sensitive_layer_index=1,sensitive_volume="GreatBoxBoard1Component",
                  suffix_adj_run_tags=["run5e5evts_splitf1_","run5e5evts_splitf10_"],
                  suffix_fwd_run="",labels=["fsplit 1","fsplit 10"],styles=["k","r"])


check_results_for_several_sims(["test_Al_1.00mmAlEq_0.01mm_cut",
                               "test_Ta_10.00mmAlEq_0.01mm_cut"],
                            ["Al 1mm","Ta 10 mm (eqAl)"],
              show_plots=True,
              save_figure=True,
              sensitive_layer_index=1,fmts=["k-","k--","r-","r--","b-","b--"])
"""





