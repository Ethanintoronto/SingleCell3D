import numpy as np
import matplotlib.pyplot as plt 
import os
from datetime import datetime
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial'] # Choose preferred fonts

dir  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
date = datetime.now().strftime("%Y-%m-%d")
#date = "2025-03-05"


period = "500"
gamma = "2"
timestep = "0p001"
id = "000"
base = "_gamma_"+gamma+"_T_"+period+"_timestep_"+timestep+"_"+id
run = "\\data\\"+date+"\\"+date+base+"\\"
file = "Single_cell_"+"Area"+base+".txt"
AreaFile = dir + run + file
VolumeFile = dir + run + "Single_cell_"+"Volume"+base+".txt"
AreaData = np.loadtxt(AreaFile, delimiter = ",")
VolumeData = np.loadtxt(VolumeFile, delimiter = ",")

area_time = AreaData[:,0]*10
area = AreaData[:,1]

volume_time = VolumeData[:,0]*10
volume = VolumeData[:,1] 
# Create the main figure and axis
fig, ax = plt.subplots()
ax.plot(volume_time, volume, color='blue')
ax.set_xlabel(r'Time ($\tau$)', fontsize = 14)
ax.set_ylabel('Volume', fontsize = 14)
ax.set_title('Uniform Facial Stiffness', fontsize = 18)

# Create the inset axis
ax_inset = inset_axes(ax, width="40%", height="30%", loc='upper right')
ax_inset.plot(area_time, area, color='red')
ax_inset.tick_params(labelsize=8)
ax_inset.set_xlabel(r'Time ($\tau$)', fontsize=12)
ax_inset.set_ylabel('Area', fontsize=12)
plt.savefig("plots\\"+ date+"\\"+date+base+"\\"+date+base+"inlaidVA")

plt.show()


