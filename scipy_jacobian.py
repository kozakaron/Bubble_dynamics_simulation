"""
This is a modified and JIT-ted version of the internal numerical Jacobian approximation used by scipy.integrate.solve_ivp, when no analytical Jacobian is provided.
The original code can be found at `scipy/integrate/_ivp/common.py` in the SciPy source code.
"""


import numpy as np
from numba import njit
from numba.types import Tuple, float64, int64   # JIT types

# To avoid unnecessary recompilation with numba JIT
fun_id = 0
num_jac_wrapped = None

EPS = np.finfo(float).eps
NUM_JAC_DIFF_REJECT = EPS ** 0.875
NUM_JAC_DIFF_SMALL = EPS ** 0.75
NUM_JAC_DIFF_BIG = EPS ** 0.25
NUM_JAC_MIN_FACTOR = 1e3 * EPS
NUM_JAC_MAX_FACTOR = 1e10

def wrap_num_jac(fun):
    def num_jac(t, y, threshold, factor, P_amb, alfa_M, T_inf, surfactant, P_v, mu_L, rho_L, c_L, ex_args, extra_dims):
        neval = 0
        ndim = len(y)
        f_0 = fun(t, y, P_amb, alfa_M, T_inf, surfactant, P_v, mu_L, rho_L, c_L, ex_args, extra_dims)
        neval += 1

    # Calculate step size (h)
        f_sign = 2.0 * (f_0 >= 0.0).astype(np.float64) - 1.0
        y_scale = f_sign * np.maximum(threshold, np.abs(y))
        h = (y + factor * y_scale) - y  # y+(...)-y to ensure numerical stability
        # make sure that h != 0
        for i in np.nonzero(h == 0)[0]:
            while h[i] == 0 and factor[i] < NUM_JAC_MAX_FACTOR:
                factor[i] *= 10.0
                h[i] = (y[i] + factor[i] * y_scale[i]) - y[i]

    # Calculate numerical Jacobian
        f_new = np.empty((ndim, ndim), dtype=np.float64)
        for i in range(ndim):
            y_temp = y.copy()
            y_temp[i] += h[i]
            f_new[:, i] = fun(t, y_temp, P_amb, alfa_M, T_inf, surfactant, P_v, mu_L, rho_L, c_L, ex_args, extra_dims)
            neval += 1
        diff = f_new - f_0[:, None]
    # Check if the numerical Jacobian is acceptable
        max_ind = np.argmax(np.abs(diff), axis=0)
        max_diff = np.empty((ndim), dtype=np.float64)
        for i in range(ndim):
            max_ind_i = max_ind[i]
            max_diff[i] = np.abs(diff[max_ind_i, i])
        
        scale = np.empty((ndim), dtype=np.float64)
        for i in range(ndim):
            max_ind_i = max_ind[i]
            scale[i] = np.maximum(np.abs(f_0[max_ind_i]), np.abs(f_new[max_ind_i, i]))
        diff_too_small = max_diff < NUM_JAC_DIFF_REJECT * scale

    # Recalculate the numerical Jacobian where unacceptable
        if np.any(diff_too_small):
            ind = np.nonzero(diff_too_small)[0]
            new_factor = 10.0 * factor[ind]
            h_new = (y[ind] + new_factor * y_scale[ind]) - y[ind]
            f_new = np.empty((ndim, len(ind)), dtype=np.float64)
            for i, j in enumerate(ind):
                y_temp = y.copy()
                y_temp[j] += h_new[i]
                f_new[:, i] = fun(t, y_temp, P_amb, alfa_M, T_inf, surfactant, P_v, mu_L, rho_L, c_L, ex_args, extra_dims)
                neval += 1
            diff_new = f_new - f_0[:, None]
    # Check if the numerical Jacobian is acceptable
            max_ind = np.argmax(np.abs(diff_new), axis=0)
            r = np.arange(len(ind))
            max_diff_new = np.empty((len(ind)), dtype=np.float64)
            for i in range(len(ind)):
                max_ind_i = max_ind[i]
                max_diff_new[i] = np.abs(diff_new[max_ind_i, i])
            scale_new = np.empty((len(ind)), dtype=np.float64)
            for i in range(len(ind)):
                max_ind_i = max_ind[i]
                scale_new[i] = np.maximum(np.abs(f_0[max_ind_i]), np.abs(f_new[max_ind_i, r[i]]))
            update = np.empty((len(ind)), dtype=np.bool_)
            for i in range(len(ind)):
                update[i] = max_diff[ind[i]] * scale_new[i] < max_diff_new[i] * scale[ind[i]]

    # Update the numerical Jacobian if needed
            if np.any(update):
                update = np.nonzero(update)[0]
                update_ind = ind[update]
                for i in range(len(update)):
                    update_ind_i = update_ind[i]
                    factor[update_ind_i] = new_factor[update[i]]
                    h[update_ind_i] = h_new[update[i]]
                    diff[:, update_ind_i] = diff_new[:, update[i]]
                    scale[update_ind_i] = scale_new[update[i]]
                    max_diff[update_ind_i] = max_diff_new[update[i]]

    # update factor
        diff /= h
        factor[max_diff < NUM_JAC_DIFF_SMALL * scale] *= 10.0
        factor[max_diff > NUM_JAC_DIFF_BIG * scale] *= 0.1
        factor = np.maximum(factor, NUM_JAC_MIN_FACTOR)
        factor = np.minimum(factor, NUM_JAC_MAX_FACTOR)
        
        return diff, factor, neval

    try:                       # ret:  diff, factor, neval
        num_jac_jit = njit(    # args: t, y, threshold, factor, P_amb, alfa_M, T_inf, surfactant, P_v, mu_L, rho_L, c_L, ex_args, extra_dims
            Tuple((float64[:, :], float64[:], int64))(float64, float64[:], float64, float64[:], float64, float64, float64, float64, float64, float64, float64, float64, float64[:], int64)
            )(
                num_jac
            )
    except Exception as error:
        print('ODE function cannot be numba JIT-ted, thus the numerical Jacobian cannot be JIT-ted.')
        raise error

    return num_jac_jit

class Jacobian:
    def __init__(self,
                 fun,
                 ndim: int,
                 threshold: float = 1e-10
    ):
        self.fun = fun
        self.ndim = ndim
        self.threshold = threshold
        self.njac = 0
        self.nfev_jac = 0
        self.factor = np.full(self.ndim, EPS ** 0.5)
        self.args = self.fun.args

        global fun_id
        global num_jac_wrapped
        if fun.id != fun_id:
            num_jac_wrapped = wrap_num_jac(fun.raw_fun)
            fun_id = fun.id

    def __call__(self, t: float, y: np.ndarray, f_0: np.ndarray=None):
        self.njac += 1
        J, factor, neval = num_jac_wrapped(t, y, self.threshold, self.factor, *self.args)
        self.factor = factor
        self.fun.add_to_neval(neval)
        self.nfev_jac += neval

        return J