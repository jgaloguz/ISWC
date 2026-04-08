# Import libraries
import matplotlib.pyplot as plt
import numpy as np

# Load data from file
spectrum1 = np.loadtxt("intensity_1au.dat")

# Prepare arrays
Energy1 = spectrum1[:,0]
ModulatedSpectrum1 = spectrum1[:,1]
UnmodulatedSpectrum1 = spectrum1[:,2]

# Plot data
ax = plt.figure().add_subplot()
ax.loglog(Energy1, ModulatedSpectrum1, label='Modulated')
ax.loglog(Energy1, UnmodulatedSpectrum1, label='Unmodulated')

# Set axes labels and legend
ax.set_xlabel("Energy")
ax.set_ylabel("Spectrum")
# ax.set_aspect('equal')
ax.legend()

# Show plot
plt.show()