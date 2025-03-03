import numpy as np
import matplotlib.pyplot as plt 
import os
from datetime import datetime

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial'] # Choose preferred fonts

dir  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
date = datetime.now().strftime("%Y-%m-%d")
period = "500"
gamma = "2"
timestep = "0p001"
id = "003"
base = "_gamma_"+gamma+"_T_"+period+"_timestep_"+timestep+"_"+id
run = "\\data\\"+date+"\\"+date+base+"\\"
file = "Single_cell_"+"MaxForce"+base+".txt"
file = dir + run + file
data = np.loadtxt(file, delimiter = ",")[1:,:]
time = data[:,0]

currtime = time[0]
forceDiffs = []
forceDiff = []
force = []
forces = []
uniq_times = []
for i in range(len(time)):
    if i == 0:
        force.append(data[i,1])
        continue
    if time[i] == currtime:
        force.append(data[i,1]/force[0])
    else:
        uniq_times.append(currtime)
        currtime = time[i]
        force[0] = 1
        forces.append(force)
        force = []
        force.append(data[i,1])

#Plot the number of midsteps vs the force with each uniq as a seperate line
num_midsteps = 10000
log_midsteps = 1000

#units of tau
uniq_times = np.array(uniq_times)*100
midsteps = np.arange(0, num_midsteps, log_midsteps)
plt.figure()
plt.xlabel("Number of Midsteps", fontsize = "large")
plt.ylabel("Max Force (Normalized to first midstep)",  fontsize = "large")
plt.title("Midsteps Required to Relax System")
for i in range(len(forces)):
    plt.plot(midsteps, forces[i], label = "time = " + str(uniq_times[i])+ r'$\tau$')
if not os.path.exists("plots\\"+ date+"\\"+date+base+"\\"):
    os.makedirs("plots\\"+ date+"\\"+date+base)
plt.legend()
plt.savefig("plots\\"+ date+"\\"+date+base+"\\"+"MaxForce_midsteps_"+id)
plt.show()

