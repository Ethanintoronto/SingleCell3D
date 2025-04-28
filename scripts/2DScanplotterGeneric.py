import matplotlib.pyplot as plt
import numpy as np
import os
import glob
from datetime import datetime

def convert_decimal(decimal_str):
    return float(decimal_str.replace('p', '.'))

# Set base directory and date
date = "2025-03-19"
data_dir = os.path.join("data", date)

# Define the range of IDs (adjust these values based on how IDs are structured in the directory names)
# Define the start and end of the ID range (e.g., from "000" to "050")
id_start = 1  # Starting ID (inclusive)
id_end = 400  # Ending ID (inclusive)
period = 500

id_range = range(id_start, id_end + 1)

# Define the range of x and y variables
x_var = r'$K_s$'
y_var = r'$K_{sTrailing}$'

x_start = 0.5    # Starting value for x
num_x  = 20     # Ending value for x
x_increment = 0.5  # Increment value for x
x_stop = x_start + x_increment * (num_x - 1)

y_increment = 0.5  # Increment value for y
num_y = 20
y_start = 0.5
y_stop = y_start + y_increment * (num_y - 1)


# Lists to store extracted values
x_values = []
y_values = []
displacements = []
for x in np.linspace(x_start, x_stop, num_x):

    for y in np.linspace(y_start, y_stop, num_y):

        # Use the looped values of x and y
        x_values.append(x)
        y_values.append(y)

# Find all relevant subdirectories
subdirs = glob.glob(os.path.join(data_dir, f"{date}_gamma_*_T_*_timestep_*_*"))
# Sort subdirectories by ID (assuming the ID is the last part)
subdirs.sort(key=lambda subdir: int(subdir.split("_")[-1]))  # Sort by the second-to-last part (ID)

for subdir in subdirs:

    # Skip old directories that do not contain "_T_"
    if "_T_" not in subdir:
        continue
        # Extract gamma, period, and ID from directory name
    parts = os.path.basename(subdir).split("_")
    
    try:
        id_index = parts.index("0p001") + 1  # Assuming the ID is right after "timestep"
        id_str = parts[id_index]
    
        id_val = int(id_str)  # Assuming ID is an integer

        # Check if ID is within the specified range
        if id_val not in id_range:
            continue
        
    except (ValueError, IndexError) as e:
        continue  # Skip if parsing fails
    # Find the data file in the subdirectory
    file_pattern = os.path.join(subdir, "Single_cell_Centroid_gamma_*_T_*_timestep_*_*.txt")
    files = glob.glob(file_pattern)
    if not files:
        continue  # Skip if no matching file found
    file = files[0]

    try:
        # Load data and compute displacement
        data = np.loadtxt(file, delimiter=",")
        # Get displacements starting after 4T:
        data = data[data[:, 0] > (4 * period / 1000)]
        start_position = data[0, 1:4]  # Initial (x, y, z)
        end_position = data[-1, 1:4]   # Final (x, y, z)
        displacement = np.linalg.norm(end_position - start_position)
        # Scale displacement to cell units
        # l^2 = 50um

        print(f"Computed displacement: {displacement}")
    except Exception as e:
        print(f"Error loading data from {file}: {e}")
        continue  # Skip if loading fails

    # Store displacement for the current x, y pair
    displacements.append(displacement)

            
# Convert lists to numpy arrays
if not x_values or not y_values:
    print("No valid data found. Exiting.")
    exit()

x_values = np.array(x_values)
y_values = np.array(y_values)
displacements = np.array(displacements) * 7  # Scale displacement to cell units

# Create grid for plotting
unique_x = np.unique(x_values)
unique_y = np.unique(y_values)

disp_matrix = np.zeros((len(unique_y), len(unique_x)))
xx = np.zeros((len(unique_y), len(unique_x)))
yy = np.zeros((len(unique_y), len(unique_x)))

for i, x in enumerate(unique_x):
    for j, y in enumerate(unique_y):
        mask = (x_values == x) & (y_values == y)
        if np.any(mask):
            disp_matrix[j, i] = displacements[mask][0]
            xx[j,i] = x_values[mask][0]
            yy[j,i] = y_values[mask][0]/x_values[mask][0] 
        else:
            disp_matrix[j, i] = np.nan  # Handle missing data


# Create plot
plot_dir = os.path.join("plots", date)
fig, ax = plt.subplots()
c = ax.imshow(disp_matrix, aspect='auto', origin='lower', cmap='jet',
              extent=[min(unique_x), max(unique_x), min(unique_y), max(unique_y)])
ax.set_xlabel(x_var, fontsize=18)
ax.set_ylabel(y_var, fontsize=18)
#ax.set_title("2D Scan Displacement")
cbar = fig.colorbar(c)
cbar.set_label(r'Displacement ($\mu$m)', fontsize =18)
fig.savefig(plot_dir+"\\"+date+"_"+str(id_start)+"_"+str(id_end)+"_"+y_var+"_"+x_var+"_imshow")

fig3, ax3 = plt.subplots()
scatter = ax3.scatter(y_values/x_values,displacements,c = x_values,cmap = "jet")
cbar3 = fig3.colorbar(scatter)
cbar3.set_label("Ks", fontsize=18)
ax3.set_ylabel(r'Displacement ($\mu$m)', fontsize=18)
ax3.set_xlabel(y_var+"/"+x_var, fontsize=18)
fig3.savefig(plot_dir+"\\"+date+"_"+str(id_start)+"_"+str(id_end)+"_"+y_var+"_"+x_var+"_scatter")

# Add a new plot where ks/ks_trailing is plotted vs displacement as a line for each ks
fig4, ax4 = plt.subplots(figsize = (6.4, 5.2))
unique_xs = np.unique(x_values)
for ks in unique_xs:
    if int(ks)!=ks:
        continue 
    mask = (x_values == ks) & (y_values>ks)
    sorted_indices = np.argsort(y_values[mask] / x_values[mask])
    ax4.plot((y_values[mask] / x_values[mask])[sorted_indices], displacements[mask][sorted_indices], label=f'Ks={ks}')
ax4.set_xlabel(y_var + "/" + x_var, fontsize=18)
ax4.set_ylabel(r'Displacement ($\mu$m)', fontsize=18)
ax4.legend(fontsize=10)
fig4.savefig(plot_dir+"\\"+date+"_"+str(id_start)+"_"+str(id_end)+"_"+y_var+"_"+x_var)

fig5, ax5 = plt.subplots()
unique_xs = np.unique(x_values)
for ks in unique_xs:
    mask = (x_values == ks) & (y_values>=ks)
    sorted_indices = np.argsort(y_values[mask] / x_values[mask])
    ax5.plot(y_values[mask][sorted_indices], displacements[mask][sorted_indices], label=f'Ks={ks}')
ax5.set_xlabel(y_var, fontsize=18)
ax5.set_ylabel(r'Displacement ($\mu$m)', fontsize=18)
ax5.legend(fontsize=10, loc = "upper right")
fig5.savefig(plot_dir+"\\"+date+"_"+str(id_start)+"_"+str(id_end)+"_"+y_var)
# Show the plot
plt.show()