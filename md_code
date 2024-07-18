import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

# Load PSF and DCD files
u = mda.Universe('data-p.psf', 'eq-1.dcd', format="LAMMPS")

# Select atoms named 'NIT'
ag = u.select_atoms('resname NIT')

# Calculate center of mass of selected atoms 
center_of_mass = ag.center_of_mass()

# Set parameters
max_distance = 40  
num = 0 # Number radii considered

# Initialize dictionary for counting occurrences in each bin
bins = np.arange(1, num + 1)
counts = np.zeros((len(u.trajectory), num+1))
r_count = np.zeros(num+1)
outer_volumes = np.zeros(num+1)
inner_volumes = np.zeros(num+1)
volumes = np.zeros(num+1)

# Loop through each frame 
for ts in u.trajectory[0:1]:
    for atom in ag:
        # Calculate distance from center of mass along Z-axis
        distance = np.sqrt(np.sum((atom.position - center_of_mass)**2))
        if distance <= max_distance:
            bin_index = int(distance / max_distance * num) + 1  # Calculate bin index
            if bin_index in bins:
                counts[ts.frame, bin_index] += 1

average = np.mean(counts,axis = 0)

for i in range (num +1): 
    # Calculate volumes of spherical shells (hollow spheres)
    radial_count = abs(average[i+1] - average[i])
    outer_volumes = (4/3) * np.pi * (bins[i]** 3)
    inner_volumes = (4/3) * np.pi * ((bins [i- 1]) ** 3)
    volumes[i] = outer_volumes - inner_volumes
    

# print(average[9])


# Calculate densities
densities = (average / volumes) * (((10**4)*6.023)/28)
print(densities)

# Plot the radial density profile along Z-axis
plt.figure(figsize=(8, 6))
plt.plot(bins, densities, linestyle='-')

plt.xlabel('Distance Bin')
plt.ylabel('Number Density (number of atoms per cubic Angstrom)')
plt.title('Number Density Profile along Radius')
plt.grid(True)
plt.tight_layout()
plt.show()

# Determine maximum observed density
max_density = np.max(densities)

# Set threshold density as 0.02
threshold_density = 0.05 * max_density  # 2% of maximum density

# Find distances where densities drop below threshold
for i in range (0,40): 
    if densities[i] < threshold_density and i > 25:
            print(f"The radius is approximately {i:.6f} Angstroms.")
            break;
    else:
            continue;



