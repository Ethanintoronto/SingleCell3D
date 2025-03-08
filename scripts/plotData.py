import numpy as np
import matplotlib.pyplot as plt 
import os
from datetime import datetime
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial'] # Choose preferred fonts

dir  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
date = datetime.now().strftime("%Y-%m-%d")
date = "2025-03-05"


period = "500"
gamma = "2"
timestep = "0p001"
id = "001"
base = "_gamma_"+gamma+"_T_"+period+"_timestep_"+timestep+"_"+id
run = "\\data\\"+date+"\\"+date+base+"\\"
tests =  ["Area", "Centroid", "Volume", "MaxForce", "GeoCentroid"]
plot = ["Single_cell_"+x+base+".txt" for x in tests]
filenames = [dir + run + file for file in plot] 

for index, file in enumerate(filenames):
    data = np.loadtxt(file, delimiter = ",")
    time = data[:,0]
    if tests[index] == "Centroid":
        plotLabels = [r'Centroid x ($\mu$m)', r'Centroid y ($\mu$m)',r'Centroid z ($\mu$m)']
        saveLabels = ["Centroid x", "Centroid y","Centroid z"]
        data[:,1:] = data[:,1:]*7 
    elif tests[index] == "GeoCentroid":
        plotLabels = [r'Geometric Centroid x ($\mu$m)',r'Geometric Centroid y ($\mu$m)',r'Geometric Centroid z ($\mu$m)']
        saveLabels = ["GeoCentroid x", "GeoCentroid y","GeoCentroid z"]
        data[:,1:] = data[:,1:]*7 
    elif tests[index] == "MaxForce":
        plotLabels = ["Max Force", "Max Force x","Max Force y","Max Force z"]
        saveLabels = ["Max Force", "Max Force x","Max Force y","Max Force z"]
    else:
        saveLabels= [tests[index]]
        plotLabels= [tests[index]]

    #Convert time from simulation time to units of tau
    #simulation timestep = 0.01tau 
    #time recorded = time_*timestep = time_*0.001 -> time in minutes = time recorded*100tau*0.5 tau/min

    #keeping units of tau for now tau = time*timestep*1000(time/timesteps)*1/100(tau/time):
    
    time = time*10    
    indep = data[:,1:]
    for i in range(len(indep[0])):
        plt.figure()
        plt.plot(time, indep[:,i])
        plt.xlabel(r'Time ($\tau$)', fontsize = 14)
        plt.ylabel(plotLabels[i], fontsize = 14)
        plt.title(str(saveLabels[i])+" with Stiff Boundary (10 midsteps)", fontsize = 18)
        if not os.path.exists("plots\\"+ date+"\\"+date+base+"\\"):
            os.makedirs("plots\\"+ date+"\\"+date+base)
        plt.savefig("plots\\"+ date+"\\"+date+base+"\\"+date+base+saveLabels[i])
plt.show()
