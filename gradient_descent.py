"""________________________________Libraries________________________________"""

import numpy as np   # matrices, math
import time   # runtime measurement
import random   # random number generator
import importlib   # reload changes you made
from termcolor import colored   # colored error messages

# my own file:
already_imported = 'de' in globals()
try:
    import full_bubble_model as de
except:
    try:
        import Bubble_dynamics_simulation.full_bubble_model as de
    except:
        print(colored(f'Error, \'full_bubble_model.py\' not found', 'red'))
if already_imported: importlib.reload(de)   # reload changes you made


"""________________________________Helper functions________________________________"""

def rand_point(ranges, ID=0, padding=0.001):
    '''Generates a random cpar dict within ranges. Arguments:
     * ranges: dict, ranges of the parameters ([single_value] or [min, max])
     * ID: int, ID for the point
     * padding: float, padding for the ranges (0.0 <= padding <= 0.5) to avoid points too close to the edge

    Returns:
     * point: dict, cpar with random values (warning: not a dotdict)
    '''

    if padding < 0.0 or 0.5 < padding:
        print(colored(f'Error, padding must be between 0.0 and 0.5, instead it is {padding}', 'red'))

    point = dict(ID=ID)
    for key in ranges:
        if len(ranges[key]) < 2:
            point[key] = ranges[key][0]
        elif (key == 'fractions'):
            padding_size = padding * abs(ranges[key][1][0] - ranges[key][0][0])
            fisrt_material = random.uniform(ranges[key][0][0] + padding_size, ranges[key][1][0] - padding_size) # fraction of the first material (in case of NH3 it is H2)
            point[key] = [fisrt_material, 1.0-fisrt_material]
        else:
            padding_size = padding * abs(ranges[key][1] - ranges[key][0])
            point[key] = random.uniform(ranges[key][0] + padding_size, ranges[key][1] - padding_size)

    if not 'P_v' in point:
        point['P_v'] = de.VapourPressure(point['T_inf'])  # [Pa]
    if not 'mu_L' in point:
        point['mu_L'] = de.Viscosity(point['T_inf'])  # [Pa*s]

    return point

def _in_range(point, ranges):
    '''Checks if a cpar is in ranges or not. Arguments:
     * point: dict, cpar
     * ranges: dict, ranges of the parameters ([single_value] or [min, max])
     
     Returns:
      * bool, True if the point is in ranges, False otherwise'''

    for key in ranges:
        if key in ['ID', 'gases']:
            continue
        if len(ranges[key]) < 2:
            if point[key] != ranges[key][0]:
                return False
        elif key == 'fractions': 
            if point[key][0] < ranges[key][0][0] or ranges[key][1][0] < point[key][0]:
                return False
        elif point[key] < ranges[key][0] or ranges[key][1] < point[key]:
            return False
    return True

def _squeeze_into_ranges(point, ranges, padding=0.0, verbose=False):
    '''Manages out of range points. Arguments:
     * point: dict, cpar (will be changed)
     * ranges: dict, ranges of the parameters ([single_value] or [min, max])
     * padding: float, how mouch should the point be inside the ranges (0.0 <= padding <= 0.5)
     
     Returns:
      * point: dict, cpar'''

    for key in ranges:
        if len(ranges[key]) < 2:
            continue

        if ranges[key][0] < ranges[key][1]:
            start = ranges[key][0]
            end = ranges[key][1]
        else:
            start = ranges[key][1]
            end = ranges[key][0]

        if point[key] < start:
            interval_width = abs(ranges[key][1] - ranges[key][0])
            point[key] = start + padding * interval_width
            if verbose:
                print(colored(f'\tWarning, {key} is out of range: {key}={point[key]}; min_value={ranges[key][0]}', 'yellow'))
        elif point[key] > end:
            interval_width = abs(ranges[key][1] - ranges[key][0])
            point[key] = end - padding * interval_width
            if verbose:
                print(colored(f'\tWarning, {key} is out of range: {key}={point[key]}; max_value={ranges[key][1]}', 'yellow'))

    return point

def evaluate(point, to_optimize, t_int, LSODA_timeout, Radau_timeout):
    '''Runs the simulation, and returns with data extended with output (what we want to optimize). Arguments:
     * point: dict, control parameter which we want to evaluate
     * to_optimize: str, name of the output we want to optimize (e.g. 'energy_efficiency')
     * t_int, LSODA_timeout, Radau_timeout: see de.solve() for more info
     
     Returns:
      * data: dict, simulation results from de.get_data()
      * success: bool, True if the simulation was successful, False otherwise'''
    
    # run simulation
    if not 'P_v' in point:
        point['P_v'] = de.VapourPressure(point['T_inf'])  # [Pa]
    if not 'mu_L' in point:
        point['mu_L'] = de.Viscosity(point['T_inf'])  # [Pa*s]
    cpar = de.dotdict(point)
    num_sol, error_code, elapsed_time = de.solve(cpar, t_int, LSODA_timeout, Radau_timeout)   # simulation without plotting
    data = de.get_data(cpar, num_sol, error_code, elapsed_time)   # post processing
    del num_sol
    _, success = de.get_errors(error_code)

    # return value to optimize
    output = data[to_optimize]
    if not success: # WARNING this part is designed for to_optimize='energy_efficiency' only
        output = 1.0e30
    elif output < 0.0 and data['expansion_work'] >= 0.0:
        output = 1.0e30
    point['output'] = output
    data['output'] = output
    return data, success

def evaluate_kwargs(kwargs):
    """Call evaluate() with a dictionary containing the arguments. Arguments:
    * kwargs: dict, containing the arguments for evaluate()
        
    Returns:
    * data: dict, simulation results from de.get_data()
    * point: dict, control parameter which we want to evaluate
    """

    point = kwargs['point']
    data, success = evaluate(**kwargs)
    point['success'] = success
    return [dict(data), dict(point), success]

def _norm_gradient(gradient, verbose=False):
    '''Calculates the norm of the gradient. Arguments:
     * gradient: dict, gradient of the point (will be changed)
     
     Returns:
      * dict, normed gradient vector'''
    
    # calculate length of the gradient
    gradient_len = 0.0
    for key in gradient:
        gradient_len += gradient[key]**2
    gradient_len = np.sqrt(gradient_len)

    # deal with zero length
    if gradient_len == 0.0:
        if verbose:
            print(colored(f'\tError, gradient length is zero: {gradient=}', 'red'))
        return None

    # norm gradient
    for key in gradient:
        gradient[key] /= gradient_len

    return gradient


"""________________________________Calculate_gradient________________________________"""

def _forward_difference(point, ranges, to_optimize='energy_efficiency', delta=1e-6, t_int=np.array([0.0, 1.0]), LSODA_timeout=30, Radau_timeout=300, verbose=True):
    '''
    Calculate the normed gradient of a point with forward difference. Arguments:
     * point: dict, cpar
     * ranges: dict, ranges of the parameters ([single_value] or [min, max])
     * to_optimize: str, name of the output we want to optimize (e.g. 'energy_efficiency')
     * delta: float, step size for the finite difference
     * t_int, LSODA_timeout, Radau_timeout: see de.solve() for more info
     * verbose: bool, print stuff

    Returns:
     * gradient: dict, normed gradient vector
     * datas: list of dicts, simulation results from de.get_data()
    '''
    datas = []
    gradient = dict()

    # central point
    central_point = de.copy(point)
    central_data, success = evaluate(central_point, to_optimize, t_int, LSODA_timeout=30, Radau_timeout=300)
    if not success:
        if verbose: print(colored(f'\tError, central point failed', 'red'))
        return None, None
    datas.append(central_data)
    
    # forward points
    for key in ranges:
        if len(ranges[key]) != 2:
            continue
        forward_point = de.copy(point)
        success = False
        current_delta = delta
        trial_num = 0
        while not success and trial_num < 3:
            value = point[key] + current_delta * abs(ranges[key][1] - ranges[key][0])
            forward_point[key] = value
            forward_data, success = evaluate(forward_point, to_optimize, t_int, LSODA_timeout, Radau_timeout)
            datas.append(forward_data)
            trial_num += 1
            if not success: # increase delta, if evaluation failed
                if verbose: print(colored(f'\tWarning, forward point failed, {trial_num=}; {current_delta=}', 'yellow'))
                current_delta *= 2
                continue

        # calculate finite difference (gradient[key])
        if success:
            gradient[key] = (forward_data['output'] - central_data['output']) / current_delta
        else:
            gradient[key] = 0.0
    
    return _norm_gradient(gradient), datas

def _central_difference(point, ranges, to_optimize='energy_efficiency', delta=1e-6, t_int=np.array([0.0, 1.0]), LSODA_timeout=30, Radau_timeout=300, verbose=True):
    '''Calculate the normed gradient of a point with central difference. Arguments:
     * point: dict, cpar
     * ranges: dict, ranges of the parameters ([single_value] or [min, max])
     * to_optimize: str, name of the output we want to optimize (e.g. 'energy_efficiency')
     * delta: float, step size for the finite difference
     * t_int, LSODA_timeout, Radau_timeout: see de.solve() for more info
     * verbose: bool, print stuff

    Returns:
     * gradient: dict, normed gradient vector
     * datas: list of dicts, simulation results from de.get_data()
    '''
    datas = []
    gradient = dict()

    for key in ranges:
        if len(ranges[key]) != 2:
            continue

        # forward point
        forward_point = de.copy(point)
        success_forward = False
        forward_delta = delta
        trial_num = 0
        while not success_forward and trial_num < 3:
            value = point[key] + forward_delta * abs(ranges[key][1] - ranges[key][0])
            forward_point[key] = value
            forward_data, success_forward = evaluate(forward_point, to_optimize, t_int, LSODA_timeout, Radau_timeout)
            trial_num += 1
            datas.append(forward_data)
            if not success_forward: # increase delta, if evaluation failed
                if verbose: print(colored(f'\tWarning, forward point failed, {trial_num=}; {forward_delta=}', 'yellow'))
                forward_delta *= 2
                continue

        # backward point
        backward_point = de.copy(point)
        success_backward = False
        backward_delta = delta
        trial_num = 0
        while success_forward and not success_backward and trial_num < 3:
            value = point[key] - backward_delta * abs(ranges[key][1] - ranges[key][0])
            backward_point[key] = value
            backward_data, success_backward = evaluate(backward_point, to_optimize, t_int, LSODA_timeout, Radau_timeout)
            trial_num += 1
            datas.append(backward_data)
            if not success_backward: # increase delta, if evaluation failed
                if verbose: print(colored(f'\tWarning, backward point failed, {trial_num=}; {backward_delta=}', 'yellow'))
                backward_delta *= 2
                continue

        # central difference
        if success_forward and success_backward:
            gradient[key] = (forward_data['output'] - backward_data['output']) / (forward_delta + backward_delta)
        else:
            gradient[key] = 0.0
            if verbose: print(colored(f'\tError, central difference failed', 'red'))
    
    return _norm_gradient(gradient), datas
    
"""________________________________Search________________________________"""

def search(kwargs):
    '''Call gradient_descent() with a dictionary containing the arguments.'''
    return gradient_descent(**kwargs)

def gradient_descent(ranges, path, to_optimize, start_point, step_limit=100, first_step=0.01, min_step=1e-4, delta=1e-6, verbose=True, t_int=np.array([0.0, 1.0]), LSODA_timeout=30, Radau_timeout=300):
    '''Gradient search. Starts at start_point, always go in the direction of the local gradient. Step size decays, if the output at the next step would
    bigger then at the current step, or if too many steps were taken without decay (to avoid back &forth stepping). Search ends, if the step_size
    decays to be smaller than min_step*interval_width, or if gradient fails repeatedly.     Arguments:
     * ranges: dict, ranges of the parameters ([single_value] or [min, max])
     * path: str, save location. If None or '', datas will be returned (not recommended due to memory usage), otherwise last data will be returned
     * to_optimize: str, name of the output we want to optimize (e.g. 'energy_efficiency')
     * start_point: dict, cpar where the search starts
     * step_limit: int, maximum number of steps
     * first_step: float, first step size
     * min_step: float, minimum step size
     * delta: float, step size for the finite difference
     * t_int, LSODA_timeout, Radau_timeout: see de.solve() for more info
     * verbose: bool, print stuff

    Returns:
     * datas: list of dicts, simulation results from de.get_data() if path is None or last data if path is given
     * last_bests: list of floats, list of the best outputs in each step
     * elapsed_time: float, elapsed time in seconds
    '''

    # setup
    if path is None or path == '':
        all_datas = []
        file = None
    else:
        file = de.Make_dir(path)
        file.new_file()
    solver_kwargs = dict(t_int=t_int, LSODA_timeout=LSODA_timeout, Radau_timeout=Radau_timeout)  # arguments for de.solve()
    keys = [key for key in ranges if len(ranges[key]) == 2]  # keys of the parameters we want to optimize

    start = time.time()
    step_size = first_step
    step_num = 0
    fail_num = 0
    start_point = de.dotdict(start_point)
    change = dict()
    
    absolute_best = 1e30
    last_best = 1e30
    last_bests = []

    # learning
    while step_num < step_limit and min_step < step_size and fail_num < 3:
        # calculate gradient (try forward difference, if it fails, try central difference, if that fails too, repeat last step)
        gradient, current_datas = _forward_difference(start_point, ranges, to_optimize, delta=delta, **solver_kwargs, verbose=verbose)
        if gradient is None:
            gradient, current_datas = _central_difference(start_point, ranges, to_optimize, delta=delta, **solver_kwargs, verbose=verbose)
        if gradient is None:
            print(colored(f'\tError, gradient can not be calculated in point {start_point}', 'red'))
            if len(change) == 0:
                return None, [1.0e30], time.time()-start
            else:
                for key in keys:
                    start_point[key] += change[key]
                start_point = _squeeze_into_ranges(start_point, ranges, padding=10.0*delta)
                step_num += 1
                fail_num += 1
                continue

        # calculate last best and print stuff
        last_best = min([data['output'] for data in current_datas])
        last_bests.append(last_best)
        absolute_best = min(absolute_best, last_best)
        if verbose: 
            print(colored(f'{step_num}. step; {last_best=: .5e}; {absolute_best=: .5e}; {step_size=: .5f}', 'green'))
            point_str = ''.join([f'{key}={start_point[key]: e}; ' for key in keys if isinstance(start_point[key], float)])
            print(f'\tpoint   =({point_str})')

        # make 5 steps in the direction of the gradient
        trial_outputs = []
        modifiers = [0.25, 0.5, 1.0, 2.0, 4.0]
        for modifier in modifiers:
            trial_point = de.copy(start_point)
            for key in keys:
                trial_point[key] -= modifier * step_size * gradient[key] * abs(ranges[key][1] - ranges[key][0])
            trial_point = _squeeze_into_ranges(trial_point, ranges, padding=10.0*delta)
            trial_data, success = evaluate(trial_point, to_optimize, **solver_kwargs)
            current_datas.append(trial_data)
            if success:
                trial_outputs.append(trial_data['output'])
            else:
                trial_outputs.append(1.0e30)

        # choose the best trial point and modify step_size accordingly
        if any([output < 1.0e30 for output in trial_outputs]):  # if there is any successfully evaluated trial point
            best_modifier = modifiers[trial_outputs.index(min(trial_outputs))]
            step_size *= best_modifier
        else:
            if verbose: print(colored(f'\tError, simulation failed in all trial points, {fail_num=}. ', 'red'))
            if len(change) == 0:
                return None, [1.0e30], time.time()-start
            else:
                for key in keys:
                    start_point[key] += change[key]
                start_point = _squeeze_into_ranges(start_point, ranges, padding=10.0*delta)
                step_num += 1
                fail_num += 1
                continue

        # make the step
        for key in keys:
            change[key] = -step_size * gradient[key] * abs(ranges[key][1] - ranges[key][0])
            start_point[key] += change[key]
        start_point = _squeeze_into_ranges(start_point, ranges, padding=10.0*delta)
        if verbose:
            next_point_data, success = evaluate(start_point, to_optimize, **solver_kwargs)
            print(f'\toutput  ={next_point_data["output"]}; success={success}; modifier={best_modifier}; ')
        
        # print stuff
        if file is not None:
            for data in current_datas:
                file.write_line(data)
                del data
        else:
            all_datas.append(current_datas)
        del current_datas

        step_num += 1
        fail_num = 0
        if verbose:
            gradient_str = ''.join([f'{key}={gradient[key]: 6.4f}; ' for key in keys])
            print(f'\tgradient=({gradient_str})')
            change_str = ''.join([f'{key}={change[key]: e}; ' for key in keys])
            print(f'\tchange  =({change_str})')

    end = time.time()
    elapsed_time = end-start
    if file is not None:
        file.close()
        return dict(trial_data), last_bests, elapsed_time
    else:
        return all_datas, last_bests, elapsed_time
