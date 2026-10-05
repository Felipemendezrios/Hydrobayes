import numpy as np
from numba.pycc import CC
from numba import njit

cc = CC('run_iqr')

@cc.export('run_iqr', 'f8[:](i4, f8[::1])')
def run_iqr(half_width, data):
    """Computes a running Innerquartile Range
    The routine accepts a column vector as input.  "halfWidth" number of data
    points for computing the Innerquartile Range are selected before and
    after the target data point, but no including the target data point.
    Near the ends of the series the number of points before or after are reduced.
    Nan in the data are counted as points.  The IQR is computed on the selected
    subset of points.  The process occurs for each point in the provided column vector.
    A column vector with the computed IQR at each point is returned.

    Parameters
    ----------
    half_width: int
        Number of ensembles before and after current ensemble which are used to compute the IQR
    data: np.array(float)
        Data for which the IQR is computed
    """
    npts = len(data)
    half_width = int(half_width)

    if npts < 20:
        half_width = int(np.floor(npts / 2))

    iqr_array = []

    # Compute IQR for each point
    for n in range(npts):

        # Sample selection for 1st point
        if n == 0:
            sample = data[1:1 + half_width]

        # Sample selection a end of data set
        elif n + half_width > npts:
            sample = np.hstack((data[n - half_width - 1:n - 1], data[n:npts]))

        # Sample selection at beginning of data set
        elif half_width >= n + 1:
            sample = np.hstack((data[0:n], data[n + 1:n + half_width + 1]))

        # Sample selection in body of data set
        else:
            sample = np.hstack((data[n - half_width:n], data[n + 1:n + half_width + 1]))

        iqr_array.append(iqr(sample))

    return np.array(iqr_array)

@njit
@cc.export('iqr', 'f8(f8[::1])')
def iqr(data_1d):
    """This function computes the iqr consistent with Matlab

    Parameters
    ----------
    data: np.ndarray
        Data for which the statistic is required

    Returns
    -------
    sp_iqr: float
        Inner quartile range

    """

    # Remove nan elements
    idx = np.where(np.logical_not(np.isnan(data_1d)))[0]
    data_1d = data_1d[idx]
    if len(data_1d) < 2:
        sp_iqr = np.nan
    else:
        # Compute statistics
        q25 = compute_quantile(data_1d, 0.25)
        q75 = compute_quantile(data_1d, 0.75)
        sp_iqr = q75 - q25

    return sp_iqr

@njit
@cc.export('compute_quantile', 'f8(f8[::1], f8)')
def compute_quantile(data_1d, q):

    sorted_data = np.sort(data_1d)
    n_samples = len(sorted_data)
    sample_idx = q * (n_samples) - 0.5
    x1 = int(np.floor(sample_idx))
    x2 = int(np.ceil(sample_idx))
    if x1 != x2:
        result = (sample_idx - x1) * (sorted_data[x2] - sorted_data[x1]) + sorted_data[x1]
    else:
        result = sorted_data[x1]
    return result

if __name__ == '__main__':
    cc.compile()