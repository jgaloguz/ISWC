# Import libraries
import matplotlib.pyplot as plt
import numpy as np

# Load data from file
trajectory1 = np.loadtxt("trajectory.dat")

# Prepare arrays x, y, z
x1 = trajectory1[:,1]
y1 = trajectory1[:,2]
z1 = trajectory1[:,3]


# Plot data
ax = plt.figure().add_subplot(projection='3d')
ax.plot(x1, y1, z1, label='Trajectory')

# Set axes labels and legend
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
# ax.set_aspect('equal')
ax.legend()

# Show plot
plt.show()