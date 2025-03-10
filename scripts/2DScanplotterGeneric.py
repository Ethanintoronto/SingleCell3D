import matplotlib.pyplot as plt
import numpy as np
import os
import glob
from datetime import datetime

def convert_decimal(decimal_str):
    return float(decimal_str.replace('p', '.'))

# Set base directory and date
date = "2025-03-10"
data_dir = os.path.join("data", date)

# Define the range of IDs (adjust these values based on how IDs are structured in the directory names)
# Define the start and end of the ID range (e.g., from "000" to "050")
id_start = 100  # Starting ID (inclusive)
id_end = 199  # Ending ID (inclusive)
period = 500

id_range = range(id_start, id_end + 1)

# Define the range of x and y variables
x_var = r'$K_s$'
y_var = r'$K_{sTrailing}$'

x_start = 0.5    # Starting value for x
num_x  = 10     # Ending value for x
x_increment = 0.5  # Increment value for x
x_stop = x_start + x_increment * (num_x - 1)

y_increment = 0.5  # Increment value for y
num_y = 10


# Lists to store extracted values
x_values = []
y_values = []
displacements = []

for x in np.linspace(x_start, x_stop, num_x):
    y_start = x
    y_stop = y_start + y_increment * (num_y - 1)
    for y in np.linspace(y_start, y_stop, num_y):

        # Use the looped values of x and y
        x_values.append(x)
        y_values.append(y)

# Find all relevant subdirectories
subdirs = glob.glob(os.path.join(data_dir, f"{date}_gamma_*_T_*_timestep_*_*"))
# Sort subdirectories by ID (assuming the ID is the last part)
subdirs.sort(key=lambda subdir: int(subdir.split("_")[-1]))  # Sort by the second-to-last part (ID)

print(f"Found {len(subdirs)} subdirectories.")

for subdir in subdirs:
    print(f"Processing: {subdir}")

    # Skip old directories that do not contain "_T_"
    if "_T_" not in subdir:
        print(f"Skipping (no _T_ found): {subdir}")
        continue
        # Extract gamma, period, and ID from directory name
    parts = os.path.basename(subdir).split("_")
    
    try:
        id_index = parts.index("0p001") + 1  # Assuming the ID is right after "timestep"
        id_str = parts[id_index]
        
        print(f"Extracted: ID={id_str}")
    
        id_val = int(id_str)  # Assuming ID is an integer

        # Check if ID is within the specified range
        if id_val not in id_range:
            print(f"Skipping (ID {id_val} not in range): {subdir}")
            continue
        
    except (ValueError, IndexError) as e:
        print(f"Skipping directory (failed parsing: {e}): {subdir}")
        continue  # Skip if parsing fails
    # Find the data file in the subdirectory
    file_pattern = os.path.join(subdir, "Single_cell_Centroid_gamma_*_T_*_timestep_*_*.txt")
    files = glob.glob(file_pattern)
    if not files:
        print(f"No data file found in: {subdir}")
        continue  # Skip if no matching file found
    file = files[0]
    print(f"Using data file: {file}")

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
fig, ax = plt.subplots()
c = ax.imshow(disp_matrix, aspect='auto', origin='lower', cmap='jet',
              extent=[min(unique_x), max(unique_x), min(unique_y), max(unique_y)])
ax.set_xlabel(x_var, fontsize="large")
ax.set_ylabel(y_var, fontsize="large")
ax.set_title("2D Scan Displacement")
fig.colorbar(c, label=r'Displacement ($\mu$m)')

fig2, ax2 = plt.subplots()
c2 = ax2.contourf(xx, yy,disp_matrix,levels = 8, cmap = 'jet')
ax2.set_xlabel(x_var, fontsize="large")
ax2.set_ylabel(y_var+"/"+x_var, fontsize="large")
ax2.set_title("2D Scan Displacement")
fig2.colorbar(c2, label=r'Displacement ($\mu$m)')
# Show the plot
plt.show()
