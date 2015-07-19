import numpy
import pandas
from scipy.interpolate import interp1d

from foldkin.coop.coop_model import CoopModelFactory
from foldkin.coop.coop_model_parameter_set import CoopModelParameterSet
from foldkin.util import convert_T_to_beta


# X_SCALE_FACTOR = 2.0 # N = 8
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


def load_Kf_file():
    df = pandas.read_csv('./logKf_vs_N.csv', header=0, index_col=0)
    return df

def compute_free_energy(N, log_Kf):
    parameters = CoopModelParameterSet()
    parameters.set_parameter('N', N)
    parameters.set_parameter('log_K_ss', -1.4318)
    parameters.set_parameter('log_K_ter', 0.2922)
    parameters.set_parameter('log_K_f', log_Kf)
    parameters.set_parameter('log_k0', 5.6)
    model_factory = CoopModelFactory()
    model = model_factory.create_model(parameters)
    boltzmann_factor_array = model.compute_boltzmann_factors()
    Q = boltzmann_factor_array.sum()
    T = 300.
    beta = convert_T_to_beta(T)
    G = -(1./beta) * numpy.log(boltzmann_factor_array / Q)
    G -= G[0]
    c = model.get_C_array() * 2.
    return c, G

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
    # crater_bottom -= 1.0 # exaggerated crater, N = 8
    # crater_bottom -= 1.5 # exaggerated crater, N = 4
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

    # b2 = b * 0.98 # N = 8
    b2 = b * 1.02 # N = 4
    # b2 = b * 0.50 # exaggerated crater, N = 8
    # b2 = b * 0.42 # exaggerated crater, N = 4
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

def compute_funnel_coordinates(N, n_radii, n_angles,
                               angle_range=(-numpy.pi*3/4., numpy.pi)):
    Kf_df = load_Kf_file()
    log_Kf = Kf_df.get_value(N, 'log_Kf')
    print N, log_Kf
    c, G = compute_free_energy(N, log_Kf)
    # print c
    # print G
    rescaled_c = rescale_x_values(c)
    rescaled_G = rescale_y_values(G)
    fit_fcn = make_fit_fcn(rescaled_c, rescaled_G)
    X, Y, Z = get_xyz_coords(fit_fcn, n_radii, n_angles, angle_range)
    return X, Y, Z

def compute_outline_coordinates(N, n_radii, angle_range):
    Kf_df = load_Kf_file()
    log_Kf = Kf_df.get_value(N, 'log_Kf')
    print N, log_Kf
    c, G = compute_free_energy(N, log_Kf)
    print c
    print G
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
