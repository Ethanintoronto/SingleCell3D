import numpy as np
import matplotlib.pyplot as plt 
import os

dir  = os.path.dirname(os.path.abspath(__file__))
os.chdir(dir)
base = "_V0_2p744_A0_11p76_timestep_0p001"
run = "\\data\\Single_cell_sim"+base+"\\"
tests =  ["Area", "Centroid", "MaxForce", "Volume"]
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
        if not os.path.exists("plots\\Single_cell_sim"+base+"\\"):
            os.makedirs("plots\\Single_cell_sim"+base)
        plt.savefig("plots\\Single_cell_sim"+base+"\\"+labels[i])

plt.show()

