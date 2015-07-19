import numpy
import pandas
from scipy.interpolate import interp1d


X_SCALE_FACTOR = 2.5 # N = 4
BASIN_START = 0.0
BASIN_END = 0.2 * X_SCALE_FACTOR
RAMP_START = BASIN_END
RAMP_END = 0.7 * X_SCALE_FACTOR
RIM_START = RAMP_END
RIM_END = 0.75 * X_SCALE_FACTOR
DOWN_SLOPE_START = RIM_END
DOWN_SLOPE_END = 0.98 * X_SCALE_FACTOR
CRATER_START = DOWN_SLOPE_END
CRATER_END = 1.0 * X_SCALE_FACTOR


def rescale_x_values(c):
    c_ramp_only = c[:-1]
    c_ramp_start = c_ramp_only[0]
    c_ramp_end = c_ramp_only[-1]
    c_ramp_range = c_ramp_end - c_ramp_start
    new_ramp_range = RAMP_END - RAMP_START
    new_c = (c_ramp_only - c_ramp_start) * (new_ramp_range/c_ramp_range)
    new_c += RAMP_START
    return new_c

def rescale_y_values(G):
    G_range = G.max() - G.min()
    new_range = 1.0
    new_G = (G - G.min()) * (new_range/G_range)
    return new_G

def make_fit_fcn(ramp, G):
    crater_bottom = G[-1]
    crater_bottom -= 1.5 # exaggerated crater
    G_ramp_only = G[:-1]
    ramp_interpolate = interp1d(ramp, G_ramp_only, kind='cubic', bounds_error=False,
                           fill_value=numpy.nan)
    # down_slope_interpolate = interp1d((DOWN_SLOPE_START, DOWN_SLOPE_END), (G_ramp_only[-1], crater_bottom),
    #                                   kind='linear', bounds_error=False, fill_value=numpy.nan)
    a = G_ramp_only[-1] + 0.05
    b = (a - crater_bottom) / (DOWN_SLOPE_END - DOWN_SLOPE_START)**2
    def down_slope_fcn(x):
        y = a - b * (x - DOWN_SLOPE_START)**2
        return y

    # b2 = b * 1.02 # N = 4
    b2 = b * 0.42 # exaggerated crater
    def rim_fcn(x):
        # y = numpy.ones_like(x) * G_ramp_only[-1]
        y = a - b2 * (x - DOWN_SLOPE_START)**2
        return y

    def fit_fcn(c):
        interpolated_G = numpy.apply_along_axis(ramp_interpolate, 0, c)
        down_slope_G = numpy.apply_along_axis(down_slope_fcn, 0, c)
        rim_G = numpy.apply_along_axis(rim_fcn, 0, c)
        combined_G = numpy.zeros_like(c)
        combined_G[numpy.where((c >= BASIN_START) & (c < BASIN_END))] = G_ramp_only[0]
        combined_G[numpy.where((c >= RAMP_START) & (c < RAMP_END))] = interpolated_G[numpy.where((c >= RAMP_START) & (c < RAMP_END))]
        combined_G[numpy.where((c >= RIM_START) & (c < RIM_END))] = rim_G[numpy.where((c >= RIM_START) & (c < RIM_END))]
        combined_G[numpy.where((c >= DOWN_SLOPE_START) & (c < DOWN_SLOPE_END))] = down_slope_G[numpy.where((c >= DOWN_SLOPE_START) & (c < DOWN_SLOPE_END))]
        combined_G[numpy.where((c >= CRATER_START) & (c <= CRATER_END))] = crater_bottom
        return combined_G
    return fit_fcn

def get_xyz_coords(fit_fcn, n_radii, n_angles, angle_range):
    r = numpy.linspace(BASIN_START, CRATER_END, n_radii)
    p = numpy.linspace(angle_range[0], angle_range[1], n_angles)
    R,P = numpy.meshgrid(r, p)
    X,Y = R * numpy.cos(P), R * numpy.sin(P)
    interpolated_c_values = CRATER_END - numpy.sqrt(X**2 + Y**2)
    interpolated_c_values[numpy.where(interpolated_c_values < 1e-6)] = BASIN_START
    Z = fit_fcn(interpolated_c_values)
    return X, Y, Z

def compute_funnel_coordinates(n_radii, n_angles,
                               angle_range=(-numpy.pi*3/4., numpy.pi)):
    c = numpy.array([0., 2., 4., 6., 8.])
    G = numpy.array([0., 1.146, 2.477, 3.891, -1.735])
    rescaled_c = rescale_x_values(c)
    rescaled_G = rescale_y_values(G)
    fit_fcn = make_fit_fcn(rescaled_c, rescaled_G)
    X, Y, Z = get_xyz_coords(fit_fcn, n_radii, n_angles, angle_range)
    return X, Y, Z

def compute_outline_coordinates(n_radii, angle_range):
    c = numpy.array([0., 2., 4., 6., 8.])
    G = numpy.array([0., 1.146, 2.477, 3.891, -1.735])
    rescaled_c = rescale_x_values(c)
    rescaled_G = rescale_y_values(G)
    fit_fcn = make_fit_fcn(rescaled_c, rescaled_G)

    r = numpy.linspace(BASIN_START, CRATER_END, n_radii)
    X = numpy.r_[r[::-1] * numpy.cos(angle_range[0]), r * numpy.cos(angle_range[1])]
    Y = numpy.r_[r[::-1] * numpy.sin(angle_range[0]), r * numpy.sin(angle_range[1])]
    interpolated_c_values = CRATER_END - numpy.sqrt(X**2 + Y**2)
    interpolated_c_values[numpy.where(interpolated_c_values < 1e-6)] = BASIN_START
    Z = fit_fcn(interpolated_c_values)
    return X, Y, Z
