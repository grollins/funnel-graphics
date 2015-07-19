from mayavi import mlab
from data_overhead import compute_funnel_coordinates, compute_outline_coordinates
from numpy import linspace, sqrt, pi

# -------
# OPTIONS 
# -------
FILENAME = 'funnel_overhead.png'
IMAGE_SIZE = (1600, 900)
BACKGROUND_COLOR = (0, 0, 0)
FUNNEL_COLORMAP = 'Blues'
HIGHLIGHT_EDGE = True
OUTLINE_COLOR = (1,1,1)
OUTLINE_THICKNESS = 10
N_RADII = 400
N_ANGLES = 200
# ANGLE_RANGE = (0, pi) # diameter slice
ANGLE_RANGE = (-pi*3/4., pi) # small cutaway

#------------
# PLOT FUNNEL
#------------
X, Y, Z = compute_funnel_coordinates(N_RADII, N_ANGLES, ANGLE_RANGE)
R = sqrt(X**2 + Y**2)
mlab.figure(bgcolor=BACKGROUND_COLOR)
m = mlab.mesh(X, Y, Z, scalars=R, line_width=4, colormap=FUNNEL_COLORMAP,
              resolution=32)

if HIGHLIGHT_EDGE:
    X2, Y2, Z2 = compute_outline_coordinates(N_RADII*10, ANGLE_RANGE)
    m2 = mlab.plot3d(X2, Y2, Z2, tube_radius=None, color=OUTLINE_COLOR,
                     line_width=OUTLINE_THICKNESS)

mlab.draw()
mlab.view(azimuth=202, elevation=65, distance=6.35) # small cutaway view
# mlab.view(azimuth=270, elevation=80, distance=7.0) # diameter slice view
mlab.savefig(FILENAME, size=IMAGE_SIZE)
# mlab.show()
