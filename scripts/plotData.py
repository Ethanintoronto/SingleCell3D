import numpy as np
import matplotlib.pyplot as plt 
import os
from datetime import datetime


dir  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
date = datetime.now().strftime("%Y-%m-%d")
V0 = "1"
A0 = "6"
gamma = "2p3"
timestep = "0p001"
base = "_gamma_"+gamma+"_V0_"+V0+"_A0_"+A0+"_timestep_"+timestep
run = "\\data\\"+date+"\\"+date+base+"\\"
tests =  ["Area", "Centroid", "MaxForce", "Volume", "Energy"]
plot = ["Single_cell_"+x+base+".txt" for x in tests]
filenames = [dir + run + file for file in plot] 

for index, file in enumerate(filenames):
    if tests[index] == "Centroid":
        labels = ["Centroid_x", "Centroid_y","Centroid_z"]
    elif tests[index] == "MaxForce":
        labels = ["MaxForce", "MaxForce_x","MaxForce_y","MaxForce_z"]
    else:
        labels = [tests[index]]
    data = np.loadtxt(file, delimiter = ",")
    time = data[:,0]
    indep = data[:,1:]
    for i in range(len(indep[0])):
        plt.figure()
        plt.plot(time, indep[:,i], label = labels[i])
        plt.xlabel("Time")
        plt.ylabel(labels[i])
        plt.title(labels[i]+base)
        plt.legend()
        if not os.path.exists("plots\\"+ date+"\\"+date+base+"\\"):
            os.makedirs("plots\\"+ date+"\\"+date+base)
        plt.savefig("plots\\"+ date+"\\"+date+base+"\\"+labels[i])
plt.show()
