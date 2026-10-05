"""bottom_discharge_extrapolation
Computes the extrapolated discharge in the bottom unmeasured portion of an ADCP transect. Methods are consistent with
equations used by TRDI and SonTek.

Example
-------

from DischargeFunctions.bottom_discharge_extrapolation import

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

cc = CC('bottom_discharge_extrapolation')


# Bottom Discharge Extrapolation with Numba
# =========================================
@cc.export('extrapolate_bot', 'f8[:](f8[:, :], b1[:, :], i8, f8, i4[:], f8[:, :], f8[:, :], f8[:], f8[:], '
                              'optional(i8), optional(f8))')
def extrapolate_bot(xprod,
                    w_valid_data,
                    transect_bot_method,
                    transect_exponent,
                    in_transect_idx,
                    depth_cell_size_m,
                    depth_cell_depth_m,
                    depth_processed_m,
                    delta_t,
                    bot_method=-1,
                    exponent=0.1667):
    """Computes the extrapolated bottom discharge

    Parameters
    ----------
    xprod: np.array(float)
        Cross product computed from the cross product method
    w_valid_data: np.array(bool)
        Valid water data
    transect_bot_method: int
        Stored bottom method (power = 0, no slip = 1)
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
    bot_method: int
        Specifies method to use for top extrapolation
    exponent: float
        Exponent to use for power extrapolation

    Returns
    -------
    q_bot: np.array(float)
        Bottom extrapolated discharge for each ensemble
    """

    # Determine extrapolation methods and exponent
    if bot_method == -1:
        bot_method = transect_bot_method
        exponent = transect_exponent

    # Use only data in transect
    w_valid_data = w_valid_data[:, in_transect_idx]
    xprod = xprod[:, in_transect_idx]
    cell_size = depth_cell_size_m[:, in_transect_idx]
    cell_depth = depth_cell_depth_m[:, in_transect_idx]
    depth_ens = depth_processed_m[in_transect_idx]
    delta_t = delta_t[in_transect_idx]

    # Compute bottom variables
    bot_rng = bot_variables(xprod, w_valid_data, cell_size, cell_depth, depth_ens)

    # Compute z
    z = np.subtract(depth_ens, cell_depth)

    # Use only valid data
    valid_data = np.logical_not(np.isnan(xprod))
    for row in range(valid_data.shape[0]):
        for col in range(valid_data.shape[1]):
            if valid_data[row, col] == False:
                z[row, col] = np.nan
                cell_size[row, col] = np.nan
                cell_depth[row, col] = np.nan

    # Compute bottom discharge
    q_bot = discharge_bot(bot_method, exponent, bot_rng, xprod,
                          cell_size, cell_depth, depth_ens, delta_t, z)

    return q_bot


@njit
@cc.export('discharge_top', 'f8[:](i8, f8, f8[:], f8[:, :], f8[:, :], f8[:, :], f8[:], f8[:], f8[:, :])')
def discharge_bot(bot_method, exponent, bot_rng, component,
                  cell_size, cell_depth, depth_ens, delta_t, z):
    """Computes the bottom extrapolated value of the provided component.

    Parameters
    ----------
    bot_method: int
        Bottom extrapolation method (Power, No Slip)
    exponent: float
        Exponent for power and no slip
    bot_rng: np.array(float)
        Range from the streambed to the bottom of the bottom most cell
    component: np.array(float)
        The variable to be extrapolated
    cell_size: np.array(float)
        Array of cell sizes (n cells x n ensembles)
    cell_depth: np.array(float)
        Depth of each cell (n cells x n ensembles)
    depth_ens: np.array(float)
        Bottom depth for each ensemble
    delta_t: np.array(float)
        Duration of each ensemble computed by QComp
    z: np.array(float)
        Relative depth from the bottom to each depth cell

    Returns
    -------
    bot_value: np.array(float)
        Total for the specified component integrated over the bottom range for each ensemble
    """

    # Initialize
    coef = np.repeat(np.nan, int(component.shape[1]))

    # Bottom power extrapolation
    if bot_method == 0:
        # Compute the coefficient for each ensemble
        # Loops are used for Numba compile purposes

        # Loop through ensembles
        for col in range(component.shape[1]):
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

    # Bottom no slip extrapolation
    elif bot_method == 1:
        # Valid data in the lower 20% of the water column or
        # the last valid depth cell are used to compute the no slip power fit
        cutoff_depth = 0.8 * depth_ens

        # Loop through the ensembles
        for col in range(cell_depth.shape[1]):
            numerator = 0.0
            denominator = 0.0
            numerator_valid = False
            denominator_valid = False
            cells_below_cutoff = False
            last_cell_depth = np.nan
            last_cell_size = np.nan
            last_z = np.nan
            last_component = np.nan

            # Verify there are valid depth cutoffs
            if np.any(np.logical_not(np.isnan(cutoff_depth))):

                # Loop through depth cells
                for row in range(cell_depth.shape[0]):

                    # Identify last valid cell by end of loop
                    if np.logical_not(np.isnan(cell_depth[row, col])):
                        last_cell_depth = cell_depth[row, col]
                        last_cell_size = cell_size[row, col]
                        last_z = z[row, col]
                        last_component = component[row, col]

                        # Use all depth cells below the cutoff (1 per loop)
                        if (cell_depth[row, col] - cutoff_depth[col]) >= 0:
                            cells_below_cutoff = True

                            # Compute numerator
                            numerator_temp = component[row, col] * cell_size[row, col]
                            if np.logical_not(np.isnan(numerator_temp)):
                                numerator_valid = True
                                numerator = numerator + numerator_temp

                                # If numerator computed, compute denominator
                                denominator_temp = ((z[row, col] + 0.5 * cell_size[row, col]) ** (exponent + 1)) \
                                                   - ((z[row, col] - 0.5 * cell_size[row, col]) ** (exponent + 1))
                                if np.logical_not(np.isnan(denominator_temp)) and denominator_temp != 0:
                                    denominator_valid = True
                                    denominator = denominator + denominator_temp

                # If there are not cells below the cutoff, use the last valid depth cell
                if np.logical_not(cells_below_cutoff):
                    if np.logical_not(np.isnan(last_cell_depth)):
                        # Compute numerator
                        numerator_temp = last_component * last_cell_size
                        if np.logical_not(np.isnan(numerator_temp)):
                            numerator_valid = True
                            numerator = numerator + numerator_temp

                            # If numerator computed, compute denominator
                            denominator_temp = ((last_z + 0.5 * last_cell_size) ** (exponent + 1)) \
                                               - ((last_z - 0.5 * last_cell_size) ** (exponent + 1))
                            if np.logical_not(np.isnan(denominator_temp)) and denominator_temp != 0:
                                denominator_valid = True
                                denominator = denominator + denominator_temp

                # If both numerator and denominator are valid compute the coefficient
                if numerator_valid and denominator_valid:
                    coef[col] = (numerator * (1 + exponent)) / denominator

    # Compute the bottom discharge of each profile
    bot_value = delta_t * (coef / (exponent + 1)) * (bot_rng**(exponent + 1))

    return bot_value


@njit
@cc.export('top_variables', 'f8[:](f8[:, :], b1[:, :], f8[:, :], f8[:, :], f8[:])')
def bot_variables(x_prod, w_valid_data, cell_size, cell_depth, depth_ens):
    """Computes the index to the bottom most valid cell in each ensemble and the range from
    the bottom to the bottom of the bottom most cell.

    Parameters
    ----------
    x_prod: np.array(float)
        Cross product computed from the cross product method
    w_valid_data: np.array(bool)
        Valid water data
    cell_size: np.array(float)
        Size of each depth cell in m
    cell_depth: np.array(float)
        Depth of each depth cell in m
    depth_ens: np.array(float)
        Processed depth for each ensemble

    Returns
    -------
    idx_bot: np.array(int)
        Index to the bottom most valid depth cell in each ensemble
    bot_rng: np.array(float)
        Range from the streambed to the bottom of the bottom most cell
    """

    # Identify valid data
    valid_data1 = np.copy(w_valid_data)
    valid_data2 = np.logical_not(np.isnan(x_prod))
    valid_data = np.logical_and(valid_data1, valid_data2)

    # Preallocate variables
    n_ensembles = int(valid_data.shape[1])
    bot_rng = np.repeat(np.nan, n_ensembles)

    # Loop through each ensemble
    for n in range(n_ensembles):

        # Identifying bottom most valid cell
        idx_temp = np.where(np.logical_not(np.isnan(x_prod[:, n])))[0]
        if len(idx_temp) > 0:
            idx_bot = idx_temp[-1]
            # Compute bottom range
            bot_rng[n] = depth_ens[n] - cell_depth[idx_bot, n] - 0.5 * cell_size[idx_bot, n]
        else:
            bot_rng[n] = 0

    return bot_rng


if __name__ == '__main__':
    # Used to compile code
    cc.compile()