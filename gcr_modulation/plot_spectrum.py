# Import libraries
import matplotlib.pyplot as plt
import numpy as np

# Load data from file
spectrum = np.loadtxt("spectrum.dat")

# Prepare arrays
Energy = spectrum[:,0]
ModulatedSpectrum = spectrum[:,1]
UnmodulatedSpectrum = spectrum[:,2]


# Plot data
ax = plt.figure().add_subplot()
ax.loglog(Energy, ModulatedSpectrum, label='Modulated')
ax.loglog(Energy, UnmodulatedSpectrum, label='Unmodulated')

# Set axes labels and legend
ax.set_xlabel("Energy")
ax.set_ylabel("Spectrum")
# ax.set_aspect('equal')
ax.legend()

# Show plot
plt.show()