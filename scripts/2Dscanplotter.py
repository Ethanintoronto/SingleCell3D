import matplotlib.pyplot as plt
import numpy as np
import os
import glob
from datetime import datetime
# Set base directory and date
date = "2025-03-03"
data_dir = os.path.join("data", date)

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
    
    # Extract gamma and period from directory name
    parts = os.path.basename(subdir).split("_")
    print(f"Split parts: {parts}")
    
    try:
        gamma_index = parts.index("gamma") + 1
        period_index = parts.index("T") + 1
        gamma_str = parts[gamma_index]
        period_str = parts[period_index]
        print(f"Extracted: Gamma={gamma_str}, Period={period_str}")
        gamma = float(gamma_str)
        period = float(period_str)
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
        start_position = data[0, 1:4]  # Initial (x, y, z)
        end_position = data[-1, 1:4]   # Final (x, y, z)
        displacement = np.linalg.norm(end_position - start_position)
        #scale displacement to cell units 
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

periods = np.array(periods)/100
gammas = np.array(gammas)
displacements = np.array(displacements)*7

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
#ax.scatter(periods,gammas,c = displacements, s= 200, cmap = 'jet')
ax.set_xlabel(r'T ($\tau$)', fontsize = "large")
ax.set_ylabel(r'$\gamma_0$',fontsize = "large" )
ax.set_title("Cell Displacement")
fig.colorbar(c, label=r'Displacement ($\mu$m)')
plt.show()