# Plot 1D animation

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
#plt.rcParams['animation.ffmpeg_path'] = '/usr/local/bin/ffmpeg'

# Animation parameters
datafilename = "trajectory.dat"                 # name of data file
videofilename = "animation.mp4"                 # name of video file
Tr = 20.0                                       # real animation time (in seconds)
stride = 100                                    # number of points per frame

# Import solution from file
datafile = open(datafilename, 'r')              # Open file in read mode
Nt = len(datafile.readlines())                  # Count number of points
datafile.close()                                # Close file

X = np.zeros(Nt)                                # x coordinate
Y = np.zeros(Nt)                                # y coordinate
Z = np.zeros(Nt)                                # y coordinate
datafile = open(datafilename, 'r')              # Open file in read mode
for i in range(Nt):                             # import ith frame
   data = datafile.readline().split()
   X[i] = float(data[1])
   Y[i] = float(data[2])
   Z[i] = float(data[3])
datafile.close()                                # Close file

# First set up the figure, the axis, and the plot element(s) we want to animate
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

curve, = ax.plot(X, Y, Z, linewidth=1.0, label='Trajectory')
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.set_xlim(-10.0,10.0)
ax.set_ylim(-3.0,10.0)
ax.set_zlim(-4.0,4.0)
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_zticklabels([])
ax.view_init(elev=70, azim=250)
ax.set_aspect('equal')
ax.legend()

# initialization function: plot the background of each frame
def init():
   curve.set_data([], [])
   curve.set_3d_properties([])
   return curve,

# animation function:  This is called sequentially
def animate(t):
   global stride, X, Y, Z
   print(t)
   curve.set_data(X[:t*stride],Y[:t*stride])
   curve.set_3d_properties(Z[:t*stride])
   return curve,

# call the animator
# interval (int) controls the delay b/w frames in ms
# repeat_delay (int) controls the delay b/w consecutive animation displays in ms
# repeat (bool) controls whether the animation display repeats or not
# blit (bool) controls whether blitting is used to optimze drawing or not
anim = animation.FuncAnimation(fig, animate, frames=Nt//stride, init_func=init,
                               interval=200, repeat_delay=0, repeat=True, blit=False)

# save the animation as an mp4
# this requires ffmpeg or mencoder to be installed
mp4writer = animation.FFMpegWriter(fps = Nt//stride//Tr)
anim.save(videofilename, writer=mp4writer, dpi=300)

# plt.show()
