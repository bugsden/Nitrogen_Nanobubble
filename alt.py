import MDAnalysis as mda    #import MDAnalysis imported as mda
import numpy as np
import matplotlib.pyplot as plt

# Load PSF and DCD files
#create a Universe
u = mda.Universe('data-p.psf', 'eq-1.dcd', format="LAMMPS")

# Select atoms named 'NIT'
ag = u.select_atoms('resname NIT')

# Calculate center of mass of selected atoms 
center_of_mass = ag.center_of_mass()

# Set parameters
max_distance = 40  
bins = 100
distance = np.zeros(bins)
counts = np.zeros(bins)

# Loop through each frame 
for ts in u.trajectory:
    for atom in ag:
        # Check if atom's Z-coordinate is greater than center of mass's Z-coordinate
        if atom.position[2] > center_of_mass[2]:
            # Calculate distance from center of mass along Z-axis
            distance = atom.position[2] - center_of_mass[2]
            bin_edges = np.histogram_bin_edges(distance, bins, range=[0,max_distance])
            bin_centers = (bin_edges[1:] + bin_edges[:-1]) / 2
            bin_index = np.searchsorted(bin_edges, distance, side='left') - 1 #searching to find the index where distance fits in the bin_edge
            if 0 <= bin_index:
                counts[bin_index] += 1
 
        # Calculate density in each bin
        bin_volumes = (4/3)*np.pi * ((bin_edges[1:] - bin_edges[:-1])**3)
        densities = counts/bin_volumes

# Plot the radial density profile along Z-axis
plt.figure(figsize=(8, 6))
plt.plot(bin_centers, densities, linestyle='-')
plt.xlabel('Distance (Angstroms)')
plt.ylabel('Number Density (number of atoms per cubic Angstrom)')
plt.title('Number Density Profile along Radius')
plt.grid(True)
plt.tight_layout()
plt.show()

# Determine maximum observed density
original_density = np.max(densities)

# Set threshold density as 0.02
threshold_density = 0.02 * original_density  # 2% of original density

# Find the Z coordinate where density drops below threshold
z_coordinate = None
for i, density in enumerate(densities):
    if density < threshold_density:
        z_coordinate = bin_centers[i]
        break

if z_coordinate is not None:
    print(f"The radius is approximately {z_coordinate:.6f} Angstroms.")
else:
    print("No Z coordinate found where density drops below 2% of original density.")
