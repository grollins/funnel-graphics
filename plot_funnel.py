from mayavi import mlab
from data import compute_funnel_coordinates, compute_outline_coordinates
from numpy import linspace, sqrt, pi

N = 4
n_radii = 200
n_angles = 100
# angle_range = (-pi*3/4., pi) # small cutaway
angle_range = (0, pi) # diameter slice
X, Y, Z = compute_funnel_coordinates(N, n_radii, n_angles, angle_range)
R = sqrt(X**2 + Y**2)
# mlab.figure(bgcolor=(1, 1, 1))
m = mlab.mesh(X, Y, Z, scalars=R, line_width=4, colormap='Blues', resolution=32)
X2, Y2, Z2 = compute_outline_coordinates(N, n_radii*4, angle_range)
m2 = mlab.plot3d(X2, Y2, Z2, tube_radius=None, color=(0,0,0), line_width=5)

# some kind of bug here...setting y visible and z invisible gives the intended
# result of showing only x and z axis labels
# mlab.axes(x_axis_visibility=True, y_axis_visibility=True, z_axis_visibility=False,
#           xlabel="Secondary Structure Formation", zlabel="Free Energy",
#           nb_labels=0)

mlab.draw()
# mlab.view(202, 65, 5.05) # small cutaway view, N = 8
# mlab.view(202, 65, 6.35) # small cutaway view, N = 4
mlab.view(270, 80, 5.65) # diameter slice view
# mlab.savefig('geoff_funnel_c.png', size=(1680, 910))
# mlab.show()


# Retrieve the LUT of the surf object.
# lut = m.module_manager.scalar_lut_manager.lut.table.to_array()

# The lut is a 255x4 array, with the columns representing RGBA
# (red, green, blue, alpha) coded with integers going from 0 to 255.

# light_blue = (120, 162, 185, 255)
# dark_blue = (0, 80, 159, 255)
# r = linspace(dark_blue[0], light_blue[0], 256)
# g = linspace(dark_blue[1], light_blue[1], 256)
# b = linspace(dark_blue[2], light_blue[2], 256)
# lut[:,0] = r
# lut[:,1] = g
# lut[:,2] = b
# for idx in xrange(128):
#     lut[idx,:] = light_blue
# for idx in xrange(128,256):
#     lut[idx,:] = dark_blue
# print lut[128,:]
# and finally we put this LUT back in the surface object. We could have
# added any 255*4 array rather than modifying an existing LUT.
# m.module_manager.scalar_lut_manager.lut.table = lut