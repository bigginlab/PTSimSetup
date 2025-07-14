import numpy as np

# used below
def switching_function(x, d0=0.0, r0=3.0):  # as per https://doi.org/10.1073/pnas.2113141118
    """
    Continuously differentiable switching function for water coordination.
    
    Parameters
    ----------
    x : float or np.ndarray
        Distance(s) between water oxygen and principal curve bead.
    d0 : float
        Inner cutoff below which f_sw = 1.
    r0 : float
        Smoothing range after d0 (i.e., outer cutoff is d0 + r0).
    Returns
    -------
    np.ndarray
        Values of the switching function.
    """
    x = np.asarray(x)
    f = np.zeros_like(x)
    # region where f = 1
    f[x < d0] = 1.0
    # region where function decays
    mask = (x >= d0) & (x <= d0 + r0)
    z = (x[mask] - d0) / r0 - 0.5
    f[mask] = z * (2 * z**2 - 1.5) + 0.5
    return f

# calculates water connectivity cv
def calculate_pcc(cec_heavy_atoms, principal_curve):    

    # docstring goes here

    # calculate coordination numbers given the principal curve
    coordination_numbers = []      
    for point in principal_curve:
        cn = 0.0  # initialize coordination number for this bead
        for heavy_atom in cec_heavy_atoms:
            r = np.linalg.norm(heavy_atom.position - point)
            cn += switching_function(r)
        coordination_numbers.append(cn)
    
    # convert coordination numbers of occupancy values
    occupancies = np.zeros(len(coordination_numbers))
    for i in range(len(coordination_numbers)):
        occupancies[i] = 1/(1 + np.exp((-coordination_numbers[i] + 1.25)/0.3))

    # pairwise connectivity calculation (between adjacent points)
    pairwise_connectivity = np.zeros(len(coordination_numbers) - 1)
    for i in range(len(coordination_numbers)-1):
        pairwise_connectivity[i] = (occupancies[i] + occupancies[i+1]) / 2

    # calculate pcc (geometric mean of pairwise connectivities)
    pcc = 1.0
    for i in range(len(pairwise_connectivity)):
        pcc *= pairwise_connectivity[i]
    pcc = pcc**(1.0/len(pairwise_connectivity))

    return pcc

# calculates water connectivity cv
def calculate_cec_weighted_pcc(cec_heavy_atoms, principal_curve, CEC):    

    # docstring goes here

    # calculate coordination numbers given the principal curve
    coordination_numbers = []      
    for point in principal_curve:
        cn = 0.0  # initialize coordination number for this bead
        for heavy_atom in cec_heavy_atoms:
            r = np.linalg.norm(heavy_atom.position - point)
            cn += switching_function(r)
        coordination_numbers.append(cn)
    
    # convert coordination numbers of occupancy values
    occupancies = np.zeros(len(coordination_numbers))
    for i in range(len(coordination_numbers)):
        occupancies[i] = 1/(1 + np.exp((-coordination_numbers[i] + 1.25)/0.3))

    # calcualte screening exponents based on proximity to CEC
    screening_exponents = np.zeros(len(coordination_numbers))
    for i in range(len(coordination_numbers)):
        screening_exponents[i] = switching_function(np.linalg.norm(CEC - principal_curve[i]))        # = fsw(|rcec-xi|)

    # pairwise connectivity calculation (between adjacent points)
    pairwise_connectivity = np.zeros(len(coordination_numbers) - 1)
    for i in range(len(coordination_numbers)-1):
        pairwise_connectivity[i] = ((occupancies[i]**screening_exponents[i]) + (occupancies[i+1]**screening_exponents[i+1])) / 2

    # calculate pcc (geometric mean of pairwise connectivities)
    pcc = 1.0
    for i in range(len(pairwise_connectivity)):
        pcc *= pairwise_connectivity[i]
    pcc = pcc**(1.0/len(pairwise_connectivity))

    return pcc

# for optimization of principal curve fitting
import numpy as np
from scipy.optimize import minimize
from scipy.interpolate import interp1d

def fit_polynomial_to_point_cloud(point_cloud, endpoints, n_points, degree=3):
    """
    Fit a smooth polynomial curve to a 3D point cloud constrained to start and end at given endpoints.
    
    Parameters:
        point_cloud: (N, 3) array of 3D points
        endpoints: tuple of two (3,) arrays, (start, end)
        n_points: number of evenly spaced output points
        degree: polynomial degree for fitting (default: 5)

    Returns:
        (n_points, 3) array of 3D coordinates along the fitted curve
    """
    start, end = np.array(endpoints[0]), np.array(endpoints[1])
    
    # Step 1: Project points onto a line between start and end
    line_vec = end - start
    line_unit = line_vec / np.linalg.norm(line_vec)
    projections = ((point_cloud - start) @ line_unit)  # scalar projections
    sort_idx = np.argsort(projections)
    t = projections[sort_idx]  # 1D parameter along line
    t = (t - t.min()) / (t.max() - t.min())  # normalize to [0,1]

    point_cloud_sorted = point_cloud[sort_idx]

    # Step 2: Fit polynomial for each dimension with endpoint constraints
    def fit_dim(y, y0, y1):
        # Vandermonde matrix
        T = np.vander(t, degree + 1, increasing=True)

        def objective(c):
            y_fit = T @ c
            return np.sum((y_fit - y)**2)

        def constraint_start(c):
            return np.polyval(c[::-1], 0) - y0

        def constraint_end(c):
            return np.polyval(c[::-1], 1) - y1

        cons = [{'type': 'eq', 'fun': constraint_start},
                {'type': 'eq', 'fun': constraint_end}]
        
        c0 = np.polyfit(t, y, degree)
        res = minimize(objective, c0[::-1], constraints=cons)
        return res.x  # Coefficients from low to high degree

    coeffs_x = fit_dim(point_cloud_sorted[:, 0], start[0], end[0])
    coeffs_y = fit_dim(point_cloud_sorted[:, 1], start[1], end[1])
    coeffs_z = fit_dim(point_cloud_sorted[:, 2], start[2], end[2])

    # Step 3: Sample N evenly spaced points along t in [0,1]
    t_new = np.linspace(0, 1, n_points)
    T_new = np.vander(t_new, degree + 1, increasing=True)
    
    x_new = T_new @ coeffs_x
    y_new = T_new @ coeffs_y
    z_new = T_new @ coeffs_z

    return np.stack([x_new, y_new, z_new], axis=1)
