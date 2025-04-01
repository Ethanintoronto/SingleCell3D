import matplotlib.pyplot as plt
import numpy as np
import os
import glob
from datetime import datetime

def convert_decimal(decimal_str):
    return float(decimal_str.replace('p', '.'))
# Set base directory and date
date = "2025-03-22"
data_dir = os.path.join("data", date)

# Define the range of IDs (adjust these values based on how IDs are structured in the directory names)
# Define the start and end of the ID range (e.g., from "000" to "050")
id_start = 100  # Starting ID (inclusive)
id_end = 199  # Ending ID (inclusive)

id_range = range(id_start,id_end+1)

# Find all relevant subdirectories
subdirs = glob.glob(os.path.join(data_dir, f"{date}_gamma_*_T_*_timestep_*_*"))
print(f"Found {len(subdirs)} subdirectories.")

# Lists to store extracted values
periods = []
gammas = []
displacements = []

for subdir in subdirs:
    print(f"Processing: {subdir}")
    
    # Skip old directories that do not contain "_T_"
    if "_T_" not in subdir:
        print(f"Skipping (no _T_ found): {subdir}")
        continue
    
    # Extract gamma, period, and ID from directory name
    parts = os.path.basename(subdir).split("_")
    print(f"Split parts: {parts}")
    
    try:
        gamma_index = parts.index("gamma") + 1
        period_index = parts.index("T") + 1
        id_index = parts.index("0p001") + 1  # Assuming the ID is right after "timestep"
        
        gamma_str = parts[gamma_index]
        period_str = parts[period_index]
        id_str = parts[id_index]
        
        print(f"Extracted: Gamma={gamma_str}, Period={period_str}, ID={id_str}")
        
        gamma = convert_decimal(gamma_str)
        period = float(period_str)
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
        data = data[data[:, 0] > (4 * period/1000)]
        start_position = data[0, 1:4]  # Initial (x, y, z)
        end_position = data[-1, 1:4]   # Final (x, y, z)
        displacement = np.linalg.norm(end_position - start_position)
        # Scale displacement to cell units 
        # l^2 = 50um 

        print(f"Computed displacement: {displacement}")
    except Exception as e:
        print(f"Error loading data from {file}: {e}")
        continue  # Skip if loading fails
    
    # Store values 
    periods.append(period)
    gammas.append(gamma)
    displacements.append(displacement)

# Convert lists to numpy arrays
if not periods or not gammas:
    print("No valid data found. Exiting.")
    exit()

periods = np.array(periods) / 100
gammas = np.array(gammas)
displacements = np.array(displacements) * 7

# Create grid for plotting
unique_periods = np.unique(periods)
unique_gammas = np.unique(gammas)
disp_matrix = np.zeros((len(unique_gammas), len(unique_periods)))

for i, gamma in enumerate(unique_gammas):
    for j, period in enumerate(unique_periods):
        mask = (gammas == gamma) & (periods == period)
        if np.any(mask):
            disp_matrix[i, j] = displacements[mask][0]
        else:
            disp_matrix[i, j] = np.nan  # Handle missing data

fig, ax = plt.subplots()
c = ax.imshow(disp_matrix, aspect='auto', origin='lower', cmap='jet', 
              extent=[min(unique_periods), max(unique_periods), min(unique_gammas), max(unique_gammas)])
ax.set_xlabel(r'T ($\tau$)', fontsize=18)
ax.set_ylabel(r'$\gamma_0$', fontsize=18)
#ax.set_title("Cell Displacement")
cbar = fig.colorbar(c)
cbar.set_label(r'Displacement ($\mu$m)', fontsize =18)
plt.show()
