"""
This is a modified version of the scipy.integrate.solve_ivp function from the SciPy library.
The original function is located at `scipy/integrate/_ivp/ivp.py` in the SciPy source code.
"""


import inspect
import numpy as np
from scipy.integrate._ivp.bdf import BDF
from scipy.integrate._ivp.radau import Radau
from scipy.integrate._ivp.rk import RK23, RK45, DOP853
from scipy.integrate._ivp.lsoda import LSODA
from scipy.integrate._ivp.common import EPS, OdeSolution
from scipy.integrate._ivp.base import OdeSolver
import time
import traceback
import importlib
from termcolor import colored


# import scipy_jacobian
try:
    import scipy_jacobian
    importlib.reload(scipy_jacobian)
except ImportError as _error:
    try:
        from  Bubble_dynamics_simulation import scipy_jacobian
        importlib.reload(scipy_jacobian)
    except ImportError as _error:
        print(colored(f'Error, \'scipy_ivp.py\' not found', 'red'))
        raise _error
    except Exception as _error:
        print(colored(f'Error, \'scipy_ivp.py\' failed to load', 'red'))
        raise _error
except Exception as _error:
    print(colored(f'Error, \'scipy_ivp.py\' failed to load', 'red'))
    raise _error


METHODS = {'RK23': RK23,
           'RK45': RK45,
           'DOP853': DOP853,
           'Radau': Radau,
           'BDF': BDF,
           'LSODA': LSODA}

class FunctionWrapper:
    def __init__(self, fun, args):
        self.raw_fun = fun
        self.args = args
        self.neval = 0
        self.id = id(fun)

    def __call__(self, t, x):
        self.neval += 1
        return self.raw_fun(t, x, *self.args)

    def add_to_neval(self, n):
        self.neval += int(n)

class OdeResult:
    """Result of ODE integration returned by `solve_ivp`. Members are:
     * t: array, time points
     * y: array, solution values
     * nstep: int, number of successful steps
     * nfev: int, number of evaluations of the right-hand side (including nfev_jac)
     * nfev_jac: int, number of evaluations of the right-hand taken by the numerical Jacobian approximator
     * njac: int, number of evaluations of the Jacobian
     * nlu: int, number of LU decompositions
     * runtime: float, time taken to solve the ODE
     * message: str, human-readable description of the termination reason
     * success: bool, True if the solver reached the interval end or a termination event occurred
     * details: str, in case of runtime error, this field contains the traceback"""

    def __init__(self, t, y, nstep=0, nfev=0, nfev_jac=0, njac=0, nlu=0, runtime=0.0,
                message='', success=True, details='', **kwargs):

        t = np.asarray(t)
        y = np.asarray(y)
        nstep = int(nstep)
        nfev = int(nfev)
        nfev_jac = int(nfev_jac)
        njac = int(njac)
        nlu = int(nlu)
        runtime = float(runtime)
        message = str(message)
        success = bool(success)
        details = str(details)

        self.t = t
        self.y = y
        self.nstep = nstep
        self.nfev = nfev
        self.nfev_jac = nfev_jac
        self.njac = njac
        self.nlu = nlu
        self.runtime = runtime
        self.message = message
        self.success = success
        self.details = details

        # Add costum fields:
        for key, value in kwargs.items():
            setattr(self, key, value)

    def __repr__(self):
        return str(self)

    def __str__(self):
        ret =  f'{"message": >9}: {self.message}\n'
        ret += f'{"success": >9}: {self.success}\n'
        ret += f'{"runtime": >9}: {self.runtime}\n'
        ret += f'{"t": >9}: [{self.t[0]} ...  {self.t[-1]}]\n'
        table1 = str( self.y[0]).replace('\n', '\n' + '           ')
        table2 = str(self.y[-1]).replace('\n', '\n' + '           ')
        ret += f'{"y": >9}: [{table1},\n           ...\n           {table2}]\n'
        ret += f'{"nstep": >9}: {self.nstep}\n'
        ret += f'{"nfev": >9}: {self.nfev}\n'
        ret += f'{"nfev_jac": >9}: {self.nfev_jac}\n'
        ret += f'{"njac": >9}: {self.njac}\n'
        ret += f'{"nlu": >9}: {self.nlu}\n'
        if self.details != '':
            details = self.details.replace('\n', '\n' + ' ' * 9)
            ret += f'{"details": >9}: {details}\n'
        
        return ret


def solve_ivp(fun, t_span, y0, method='RK45', timeout=-1.0, args=None,
                use_builtin_jac=False, compression=0, **options):
    """Solve an initial value problem for a system of ODEs.

    This function numerically integrates a system of ordinary differential
    equations given an initial value::

        dy / dt = f(t, y)
        y(t0) = y0

    Here t is a 1-D independent variable (time), y(t) is an
    N-D vector-valued function (state), and an N-D
    vector-valued function f(t, y) determines the differential equations.
    The goal is to find y(t) approximately satisfying the differential
    equations, given an initial value y(t0)=y0.

    Some of the solvers support integration in the complex domain, but note
    that for stiff ODE solvers, the right-hand side must be
    complex-differentiable (satisfy Cauchy-Riemann equations [11]_).
    To solve a problem in the complex domain, pass y0 with a complex data type.
    Another option always available is to rewrite your problem for real and
    imaginary parts separately.

    Parameters
    ----------
    fun : callable
        Right-hand side of the system: the time derivative of the state ``y``
        at time ``t``. The calling signature is ``fun(t, y)``, where ``t`` is a
        scalar and ``y`` is an ndarray with ``len(y) = len(y0)``. Additional
        arguments need to be passed if ``args`` is used (see documentation of
        ``args`` argument). ``fun`` must return an array of the same shape as
        ``y``. See `vectorized` for more information.
    t_span : 2-member sequence
        Interval of integration (t0, tf). The solver starts with t=t0 and
        integrates until it reaches t=tf. Both t0 and tf must be floats
        or values interpretable by the float conversion function.
    y0 : array_like, shape (n,)
        Initial state. For problems in the complex domain, pass `y0` with a
        complex data type (even if the initial value is purely real).
    method : string or `OdeSolver`, optional
        Integration method to use:

            * 'RK45' (default): Explicit Runge-Kutta method of order 5(4) [1]_.
              The error is controlled assuming accuracy of the fourth-order
              method, but steps are taken using the fifth-order accurate
              formula (local extrapolation is done). A quartic interpolation
              polynomial is used for the dense output [2]_. Can be applied in
              the complex domain.
            * 'RK23': Explicit Runge-Kutta method of order 3(2) [3]_. The error
              is controlled assuming accuracy of the second-order method, but
              steps are taken using the third-order accurate formula (local
              extrapolation is done). A cubic Hermite polynomial is used for the
              dense output. Can be applied in the complex domain.
            * 'DOP853': Explicit Runge-Kutta method of order 8 [13]_.
              Python implementation of the "DOP853" algorithm originally
              written in Fortran [14]_. A 7-th order interpolation polynomial
              accurate to 7-th order is used for the dense output.
              Can be applied in the complex domain.
            * 'Radau': Implicit Runge-Kutta method of the Radau IIA family of
              order 5 [4]_. The error is controlled with a third-order accurate
              embedded formula. A cubic polynomial which satisfies the
              collocation conditions is used for the dense output.
            * 'BDF': Implicit multi-step variable-order (1 to 5) method based
              on a backward differentiation formula for the derivative
              approximation [5]_. The implementation follows the one described
              in [6]_. A quasi-constant step scheme is used and accuracy is
              enhanced using the NDF modification. Can be applied in the
              complex domain.
            * 'LSODA': Adams/BDF method with automatic stiffness detection and
              switching [7]_, [8]_. This is a wrapper of the Fortran solver
              from ODEPACK.

        Explicit Runge-Kutta methods ('RK23', 'RK45', 'DOP853') should be used
        for non-stiff problems and implicit methods ('Radau', 'BDF') for
        stiff problems [9]_. Among Runge-Kutta methods, 'DOP853' is recommended
        for solving with high precision (low values of `rtol` and `atol`).

        If not sure, first try to run 'RK45'. If it makes unusually many
        iterations, diverges, or fails, your problem is likely to be stiff and
        you should use 'Radau' or 'BDF'. 'LSODA' can also be a good universal
        choice, but it might be somewhat less convenient to work with as it
        wraps old Fortran code.

        You can also pass an arbitrary class derived from `OdeSolver` which
        implements the solver.
    args : tuple, optional
        Additional arguments to pass to the user-defined functions.  If given,
        the additional arguments are passed to all user-defined functions.
        So if, for example, `fun` has the signature ``fun(t, y, a, b, c)``,
        then `jac` (if given) and any event functions must have the same
        signature, and `args` must be a tuple of length 3.
    use_builtin_jac : bool, optional
        If True, the solver will use a built-in Jacobian approximation. Otherwise,
        the improved Jacobian will be used. Doesn't matter if `jac` is provided. Default is False.
    compression: int, optional
        The data compression of the ODE solution:
            * 0: no compression, every step is saved (default)
            * 1: compression, about 1/10 of the steps are saved, but the solution's shape is preserved visually
            * 2: only the first and last 2 steps are saved
    **options
        Options passed to a chosen solver. All options available for already
        implemented solvers are listed below.
    first_step : float or None, optional
        Initial step size. Default is `None` which means that the algorithm
        should choose.
    max_step : float, optional
        Maximum allowed step size. Default is np.inf, i.e., the step size is not
        bounded and determined solely by the solver.
    rtol, atol : float or array_like, optional
        Relative and absolute tolerances. The solver keeps the local error
        estimates less than ``atol + rtol * abs(y)``. Here `rtol` controls a
        relative accuracy (number of correct digits), while `atol` controls
        absolute accuracy (number of correct decimal places). To achieve the
        desired `rtol`, set `atol` to be smaller than the smallest value that
        can be expected from ``rtol * abs(y)`` so that `rtol` dominates the
        allowable error. If `atol` is larger than ``rtol * abs(y)`` the
        number of correct digits is not guaranteed. Conversely, to achieve the
        desired `atol` set `rtol` such that ``rtol * abs(y)`` is always smaller
        than `atol`. If components of y have different scales, it might be
        beneficial to set different `atol` values for different components by
        passing array_like with shape (n,) for `atol`. Default values are
        1e-3 for `rtol` and 1e-6 for `atol`.
    jac : array_like, sparse_matrix, callable or None, optional
        Jacobian matrix of the right-hand side of the system with respect
        to y, required by the 'Radau', 'BDF' and 'LSODA' method. The
        Jacobian matrix has shape (n, n) and its element (i, j) is equal to
        ``d f_i / d y_j``.  There are three ways to define the Jacobian:

            * If array_like or sparse_matrix, the Jacobian is assumed to
              be constant. Not supported by 'LSODA'.
            * If callable, the Jacobian is assumed to depend on both
              t and y; it will be called as ``jac(t, y)``, as necessary.
              Additional arguments have to be passed if ``args`` is
              used (see documentation of ``args`` argument).
              For 'Radau' and 'BDF' methods, the return value might be a
              sparse matrix.
            * If None (default), the Jacobian will be approximated by
              finite differences.

        It is generally recommended to provide the Jacobian rather than
        relying on a finite-difference approximation.
    jac_sparsity : array_like, sparse matrix or None, optional
        Defines a sparsity structure of the Jacobian matrix for a finite-
        difference approximation. Its shape must be (n, n). This argument
        is ignored if `jac` is not `None`. If the Jacobian has only few
        non-zero elements in *each* row, providing the sparsity structure
        will greatly speed up the computations [10]_. A zero entry means that
        a corresponding element in the Jacobian is always zero. If None
        (default), the Jacobian is assumed to be dense.
        Not supported by 'LSODA', see `lband` and `uband` instead.
    lband, uband : int or None, optional
        Parameters defining the bandwidth of the Jacobian for the 'LSODA'
        method, i.e., ``jac[i, j] != 0 only for i - lband <= j <= i + uband``.
        Default is None. Setting these requires your jac routine to return the
        Jacobian in the packed format: the returned array must have ``n``
        columns and ``uband + lband + 1`` rows in which Jacobian diagonals are
        written. Specifically ``jac_packed[uband + i - j , j] = jac[i, j]``.
        The same format is used in `scipy.linalg.solve_banded` (check for an
        illustration).  These parameters can be also used with ``jac=None`` to
        reduce the number of Jacobian elements estimated by finite differences.
    min_step : float, optional
        The minimum allowed step size for 'LSODA' method.
        By default `min_step` is zero.

    Returns
    -------
    Bunch object with the following fields defined:
    t : ndarray, shape (n_points,)
        Time points.
    y : ndarray, shape (n, n_points)
        Values of the solution at `t`.
    nstep : int
        Number of sucesfull integration steps.
    nfev : int
        Number of evaluations of the right-hand side.
    njac : int
        Number of evaluations of the Jacobian.
    nlu : int
        Number of LU decompositions.

    message : string
        Human-readable description of the termination reason.
    success : bool
        True if the solver reached the interval end or a termination event
        occurred (``status >= 0``).
    details : string
        In case of runtime error, this field contains the traceback.

    """
    start_time = time.process_time()

    if method not in METHODS and not (
            inspect.isclass(method) and issubclass(method, OdeSolver)):
        raise ValueError(f"`method` must be one of {METHODS} or OdeSolver class.")

    t0, tf = map(float, t_span)
    timeout = float(timeout)
    compression = int(compression)

    if compression not in (0, 1, 2):
        raise ValueError("`compression` must be one of 0, 1, 2.")

    # Wrap ODE function and jacobian
    if args is not None:
        try:
            _ = [*(args)]
        except TypeError as exp:
            suggestion_tuple = (
                "Supplied 'args' cannot be unpacked. Please supply `args`"
                f" as a tuple (e.g. `args=({args},)`)"
            )
            raise TypeError(suggestion_tuple) from exp

        fun = FunctionWrapper(fun, args)

        jac = options.get('jac', None)
        have_jac = True
        if callable(jac):
            have_jac = False
            options['jac'] = lambda t, x: jac(t, x, *args)
        elif not use_builtin_jac:
            have_jac = False
            atol = options.get('atol', 1e-6)
            jac = scipy_jacobian.Jacobian(fun, len(y0), atol)
            options['jac'] = jac


        if method in METHODS:
            method = METHODS[method]

    solver = method(fun, t0, y0, tf, vectorized=False, **options)

    ts = [t0]
    ys = [y0]
    limit = 0.03 * y0[0]
    last_t = t0
    last_y = y0
    last_last_y = y0
    T_max = 0.0
    collapse_time = 0.0

    status = None
    details = ''
    nstep = 0
    while status is None:
        try:
            nstep += 1
            message = solver.step()
        except Exception as error:
            status = -1
            message = 'Runtime error: ' + str(error)
            details = ''.join(traceback.format_exception(error, limit=5))
            break

        if solver.status == 'finished':
            status = 0
            message = 'The solver successfully reached the end of the integration interval.'
        elif solver.status == 'failed':
            status = -1
            break
        if timeout >= 0.0 and time.process_time() - start_time > timeout:  # Timeout
            status = -1
            message = f'Solver timed out after {timeout: .2f} secs.'
            break

        t = solver.t
        y = solver.y

        if compression != 2:    # Only R_E (y[0]) is examined
            if ((compression == 0) or
                (abs(y[0] - ys[-1][0]) > limit) or                           # Save if change is > 3%
                (nstep % 50 == 0 and t-ts[-1] > 1e-7) or                     # Save every 50th step if not too close to last
                (t-ts[-1] > 1e-3) or                                         # Save if time difference is > 1 ms
                ((y[0] - last_y[0]) * (last_y[0] - [last_last_y[0]]) < 0)):  # Save loc min/max
                    ts.append(t)
                    ys.append(y)

        # Update T_max and collapse_time (only for the bubble simulation)  
        if y[2] > T_max:
            T_max = y[2]
        if collapse_time == 0.0:
            if (y[0] - last_y[0]) > 0.0 and (last_y[0] - last_last_y[0]) < 0.0 and y[0] < y0[0]:
                collapse_time = t
        

        last_t = t
        last_last_y = last_y
        last_y = y

        # end while

    if compression == 2:
        ts.append(last_t); ts.append(t)
        ys.append(last_y); ys.append(y)
    ts = np.array(ts)
    ys = np.vstack(ys)
    end_time = time.process_time()

    return OdeResult(t=ts, y=ys,
                     nstep=nstep, nfev=fun.neval, nfev_jac=0 if have_jac else jac.nfev_jac, njac=solver.njev if have_jac else jac.njac, nlu=solver.nlu,
                     runtime=end_time - start_time, message=message, success=status >= 0, details=details, T_max=T_max, collapse_time=collapse_time)
