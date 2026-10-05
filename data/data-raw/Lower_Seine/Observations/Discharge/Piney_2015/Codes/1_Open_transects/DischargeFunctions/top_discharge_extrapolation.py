"""top_discharge_extrapolation
Computes the extrapolated discharge in the top unmeasured portion of an ADCP transect. Methods are consistent with
equations used by TRDI and SonTek.

Example
-------

from DischargeFunctions.top_discharge_extrapolation import

    trans_select = getattr(data_in.depths, data_in.depths.selected)
    num_top_method = {'Power': 0, 'Constant': 1, '3-Point': 2, None: -1}
    self.top_ens =  extrapolate_top(x_prod, data_in.w_vel.valid_data[0, :, :],
                                    num_top_method[data_in.extrap.top_method],
                                    data_in.extrap.exponent, data_in.in_transect_idx, trans_select.depth_cell_size_m,
                                    trans_select.depth_cell_depth_m, trans_select.depth_processed_m, delta_t,
                                    num_top_method[top_method], exponent)
"""

import numpy as np
from numba.pycc import CC
from numba import njit

cc = CC('top_discharge_extrapolation')


# Top Discharge Extrapolation with Numba
# ======================================
@cc.export('extrapolate_top', 'f8[:](f8[:, :], b1[:, :], i8, f8, i4[:], f8[:, :], f8[:, :], f8[:], f8[:], '
                              'optional(i8), optional(f8))')
def extrapolate_top(xprod,
                    w_valid_data,
                    transect_top_method,
                    transect_exponent,
                    in_transect_idx,
                    depth_cell_size_m,
                    depth_cell_depth_m,
                    depth_processed_m,
                    delta_t,
                    top_method=-1,
                    exponent=0.1667):
    """Computes the extrapolated top discharge.

    Parameters
    ----------
    xprod: np.array(float)
        Cross product computed from the cross product method
    w_valid_data: np.array(bool)
        Valid water data
    transect_top_method: int
        Stored top method (power = 0, constant = 1, 3-point = 2)
    transect_exponent: float
        Exponent for power fit
    in_transect_idx: np.array(int)
        Indices of ensembles in transect to be used for discharge
    depth_cell_size_m: np.array(float)
        Size of each depth cell in m
    depth_cell_depth_m: np.array(float)
        Depth of each depth cell in m
    depth_processed_m: np.array(float)
        Depth for each ensemble in m
    delta_t: np.array(float)
        Duration of each ensemble computed from QComp
    top_method: int
        Specifies method to use for top extrapolation
    exponent: float
        Exponent to use for power extrapolation

    Returns
    -------
    q_top: np.array(float)
        Top extrapolated discharge for each ensemble
    """

    if top_method == -1:
        top_method = transect_top_method
        exponent = transect_exponent

    # Compute top variables
    idx_top, idx_top3, top_rng = top_variables(xprod, w_valid_data, depth_cell_size_m, depth_cell_depth_m)
    idx_top = idx_top[in_transect_idx]
    idx_top3 = idx_top3[:, in_transect_idx]
    top_rng = top_rng[in_transect_idx]

    # Get data from transect object
    cell_size = depth_cell_size_m[:, in_transect_idx]
    cell_depth = depth_cell_depth_m[:, in_transect_idx]
    depth_ens = depth_processed_m[in_transect_idx]

    # Compute z
    z = np.subtract(depth_ens, cell_depth)

    # Use only valid data
    valid_data = np.logical_not(np.isnan(xprod[:, in_transect_idx]))
    for row in range(valid_data.shape[0]):
        for col in range(valid_data.shape[1]):
            if valid_data[row, col] == False:
                z[row, col] = np.nan
                cell_size[row, col] = np.nan
                cell_depth[row, col] = np.nan

    # Compute top discharge
    q_top = discharge_top(top_method, exponent, idx_top, idx_top3, top_rng,
                          xprod[:, in_transect_idx], cell_size, cell_depth,
                          depth_ens, delta_t, z)

    return q_top


@njit
@cc.export('discharge_top', 'f8[:](i8, f8, i4[:], i4[:, :], f8[:], f8[:, :], f8[:, :], f8[:, :], f8[:], '
                            'f8[:], f8[:, :])')
def discharge_top(top_method, exponent, idx_top, idx_top_3, top_rng, component, cell_size, cell_depth,
                  depth_ens, delta_t, z):
    """Computes the top extrapolated value of the provided component.

    Parameters
    ----------
    top_method: int
        Top extrapolation method (Power = 0, Constant = 1, 3-Point = 2)
    exponent: float
        Exponent for the power extrapolation method
    idx_top: np.array(int)
        Index to the topmost valid depth cell in each ensemble
    idx_top_3: np.array(int)
        Index to the top 3 valid depth cells in each ensemble
    top_rng: np.array(float)
        Range from the water surface to the top of the topmost cell
    component: np.array(float)
        The variable to be extrapolated (xprod, u-velocity, v-velocity)
    cell_size: np.array(float)
        Array of cell sizes (n cells x n ensembles)
    cell_depth: np.array(float)
        Depth of each cell (n cells x n ensembles)
    depth_ens: np.array(float)
        Bottom depth for each ensemble
    delta_t: np.array(float)
        Duration of each ensemble compute by QComp
    z: np.array(float)
        Relative depth from the bottom of each depth cell computed in discharge top method

    Returns
    -------
    top_value: np.array(float)
        total for the specified component integrated over the top range
    """

    # Initialize return
    top_value = np.array([0.0])

    # Top power extrapolation
    if top_method == 0:
        coef = np.repeat(np.nan, int(component.shape[1]))

        # Compute the coefficient for each ensemble
        # Loops are used for Numba compile purposes

        # Loop through ensembles
        for col in range(component.shape[1]):
            # Initialize variables
            numerator = 0.0
            numerator_valid = False
            denominator_valid = False
            denominator = 0.0

            # Loop through depth cells in an ensemble
            for row in range(component.shape[0]):

                # Compute the numerator
                numerator_temp = component[row, col] * cell_size[row, col]
                if np.logical_not(np.isnan(numerator_temp)):
                    numerator_valid = True
                    numerator = numerator + numerator_temp

                # Compute the denominator
                denominator_temp = ((z[row, col] + 0.5 * cell_size[row, col]) ** (exponent + 1)) \
                                   - ((z[row, col] - 0.5 * cell_size[row, col]) ** (exponent + 1))
                if np.logical_not(np.isnan(denominator_temp)) and denominator_temp != 0:
                    denominator_valid = True
                    denominator = denominator + denominator_temp

            # If both numerator and denominator are valid compute the coefficient
            if numerator_valid and denominator_valid:
                coef[col] = (numerator * (1 + exponent)) / denominator

        # Compute the top discharge for each ensemble
        top_value = delta_t * (coef / (exponent + 1)) * \
            (depth_ens**(exponent + 1) - (depth_ens-top_rng)**(exponent + 1))

    # Top constant extrapolation
    elif top_method == 1:
        n_ensembles = len(delta_t)
        top_value = np.repeat(np.nan, n_ensembles)
        for j in range(n_ensembles):
            if idx_top[j] >= 0:
                top_value[j] = delta_t[j] * component[idx_top[j], j] * top_rng[j]

    # Top 3-point extrapolation
    elif top_method == 2:
        # Determine number of bins available in each profile
        valid_data = np.logical_not(np.isnan(component))
        n_bins = np.sum(valid_data, axis=0)
        # Determine number of ensembles
        n_ensembles = len(delta_t)
        # Preallocate qtop vector
        top_value = np.repeat(np.nan, n_ensembles)

        # Loop through ensembles
        for j in range(n_ensembles):

            # Set default to constant
            if (n_bins[j] < 6) and (n_bins[j] > 0) and (idx_top[j] >= 0):
                top_value[j] = delta_t[j] * component[idx_top[j], j] * top_rng[j]

            # If 6 or more bins use 3-pt at top
            if n_bins[j] > 5:
                sumd = 0.0
                sumd2 = 0.0
                sumq = 0.0
                sumqd = 0.0

                # Use loop to sum data from top 3 cells
                for k in range(3):
                    if np.isnan(cell_depth[idx_top_3[k, j], j]) == False:
                        sumd = sumd + cell_depth[idx_top_3[k, j], j]
                        sumd2 = sumd2 + cell_depth[idx_top_3[k, j], j]**2
                        sumq = sumq + component[idx_top_3[k, j], j]
                        sumqd = sumqd + (component[idx_top_3[k, j], j] * cell_depth[idx_top_3[k, j], j])
                delta = 3 * sumd2 - sumd**2
                a = (3 * sumqd - sumq * sumd) / delta
                b = (sumq * sumd2 - sumqd * sumd) / delta

                # Compute discharge for 3-pt fit
                qo = (a * top_rng[j]**2) / 2 + b * top_rng[j]
                top_value[j] = delta_t[j] * qo

    return top_value


@njit
@cc.export('top_variables', '(f8[:, :], b1[:, :], f8[:, :], f8[:, :])')
def top_variables(xprod, w_valid_data, depth_cell_size_m, depth_cell_depth_m):
    """Computes the index to the top and top three valid cells in each ensemble and
    the range from the water surface to the top of the topmost cell.

    Parameters
    ----------
    xprod: np.array(float)
        Cross product computed from the cross product method
    w_valid_data: np.array(bool)
        Valid water data
    depth_cell_size_m: np.array(float)
        Size of each depth cell in m
    depth_cell_depth_m: np.array(float)
        Depth of each depth cell in m

    Returns
    -------
    idx_top: np.array(int)
        Index to the topmost valid depth cell in each ensemble
    idx_top_3: np.array(int)
        Index to the top 3 valid depth cell in each ensemble
    top_rng: np.array(float)
        Range from the water surface to the top of the topmost cell
    """

    # Get data from transect object
    valid_data1 = np.copy(w_valid_data)
    valid_data2 = np.logical_not(np.isnan(xprod))
    valid_data = np.logical_and(valid_data1, valid_data2)

    # Preallocate variables
    # NOTE: Numba does not support np.tile
    n_ensembles = int(valid_data.shape[1])
    idx_top = np.repeat(-1, int(valid_data.shape[1]))
    idx_top_3 = np.ones((3, int(valid_data.shape[1])), dtype=np.int32)
    idx_top_3[:] = int(-1)
    top_rng = np.repeat(np.nan, n_ensembles)

    # Loop through ensembles
    for n in range(n_ensembles):
        # Identify topmost 1 and 3 valid cells
        idx_temp = np.where(np.logical_not(np.isnan(xprod[:, n])))[0]
        if len(idx_temp) > 0:
            idx_top[n] = idx_temp[0]
            if len(idx_temp) > 2:
                for k in range(3):
                    idx_top_3[k, n] = idx_temp[k]
            # Compute top range
            top_rng[n] = depth_cell_depth_m[idx_top[n], n] - 0.5 * depth_cell_size_m[idx_top[n], n]
        else:
            top_rng[n] = 0
            idx_top[n] = 0

    return idx_top, idx_top_3, top_rng


if __name__ == '__main__':
    # Used to compile code
    cc.compile()
