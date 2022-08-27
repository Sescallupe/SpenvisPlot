#path declaration and import
import sys
#It is important to set the directory, where  the subdirectory python_utilities is located 
# in the python  path this is done below
sys.path.append("..")
from python_utilities.SpenvisCSVFileHandler import SpenvisCSVFileHandler
from python_utilities.Plotters import Plotter
import pylab as pl

#Create the file handler
theSpenvisFileHandler = SpenvisCSVFileHandler()
thePlotter=Plotter()




#Test with Spenvis radiation file
##################################

file_name="../csv_files/SpenvisRadFile.csv"
proton_flux,electron_flux=theSpenvisFileHandler.ReadSpenvisTrappedRadiationFile(file_name)
print proton_flux.keys()
print proton_flux["Energy"].keys()


pl.figure(figsize=(10,10))
pl.loglog(proton_flux["Energy"]['data'],proton_flux["IFlux"]['data'],label="proton")
pl.plot(electron_flux["Energy"]['data'],electron_flux["IFlux"]['data'],label="electron")
pl.legend()
pl.xlabel("Ekin [%s]" %(proton_flux["Energy"]['unit']))
pl.ylabel("IFlux [%s]" %(proton_flux["IFlux"]['unit']))
pl.savefig("../plots/SpenvisTrapRadExample.pdf")
pl.show()

#Test with GRAS files
##################################
file_name="../csv_files/tests_analysis_niel.csv"
ModuleDescription, StatDoubleTable, Histo1Ds=theSpenvisFileHandler.ReadGRASFile(file_name)
print StatDoubleTable
print Histo1Ds['niel1'].keys()
histo1D=Histo1Ds['niel1']["NIEL SPECTRUM"]
#histo1D=Histo1Ds['source1']["KINE DISTRIBUTION"]
print histo1D
thePlotter.PlotGRAS1DHisto(histo1D,xlabel=None,
                            ylabel=None,xlog=None,
                            ylog=None,new_plot=False,histo_style=True,label=None)
pl.show()



