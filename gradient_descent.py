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

def rand_point(ranges, ID=0, padding=0.1):
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

def in_range(point, ranges):
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

def squeeze_into_ranges(point, ranges, padding=0.0, verbose=False):
    '''Manages out of range points. Arguments:
     * point: dict, cpar (will be changed)
     * ranges: dict, ranges of the parameters ([single_value] or [min, max])
     * padding: float, how mouch should the point be inside the ranges (0.0 <= padding <= 0.5)
     
     Returns:
      * point: dict, cpar'''

    for key in ranges:
        if len(ranges[key]) < 2:
            continue
        if point[key] < ranges[key][0]:
            interval_width = abs(ranges[key][1] - ranges[key][0])
            point[key] = ranges[key][0] + padding * interval_width
            if verbose:
                print(colored(f'\tWarning, {key} is out of range: {key}={point[key]}; min_value={ranges[key][0]}', 'yellow'))
        elif point[key] > ranges[key][1]:
            interval_width = abs(ranges[key][1] - ranges[key][0])
            point[key] = ranges[key][1] - padding * interval_width
            if verbose:
                print(colored(f'\tWarning, {key} is out of range: {key}={point[key]}; max_value={ranges[key][1]}', 'yellow'))

    return point

def evaluate(point, to_optimize, t_int, LSODA_timeout, Radau_timeout, log10=False):
    '''Runs the simulation, and returns with data extended with output (what we want to optimize). Arguments:
     * point: dict, control parameter which we want to evaluate
     * to_optimize: str, name of the output we want to optimize (e.g. 'energy_efficiency')
     * t_int, LSODA_timeout, Radau_timeout: see de.solve() for more info
     * log10: bool, True if we want to optimize the log10 of the output, False otherwise
     
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
    _, success = de.get_errors(error_code)

    # return value to optimize
    output = data[to_optimize]
    if not success: # WARNING this part is designed for to_optimize='energy_efficiency' only
        output = 1.0e30
    elif output < 0.0 and data['expansion_work'] >= 0.0:
        output = 1.0e30
    if log10:
        output = np.log10(output)
    point['output'] = output
    data['output'] = output
    return data, success

def norm_gradient(gradient, verbose=False):
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

# TODO implement verbose
def forward_difference(point, ranges, to_optimize='energy_efficiency', delta=1e-6, log10=True, t_int=np.array([0.0, 1.0]), LSODA_timeout=30, Radau_timeout=300):
    #keys = ['R_E', 'ratio', 'P_amb', 'alfa_M', 'T_inf', 'surfactant'] + de.excitation_args
    datas = []
    gradient = dict()

    # central point
    central_point = de.copy(point)
    central_data, success = evaluate(central_point, to_optimize, t_int, LSODA_timeout=30, Radau_timeout=300, log10=log10)
    if not success:
        print(colored(f'\tError, central point failed', 'red'))
        return None, None
    datas.append(central_data)
    
    # forward points
    for key in ranges:
        if len(ranges[key]) < 2:
            continue
        forward_point = de.copy(point)
        success = False
        current_delta = delta
        trial_num = 0
        while not success and trial_num < 3:
            value = point[key] + current_delta * abs(ranges[key][1] - ranges[key][0])
            forward_point[key] = value
            forward_data, success = evaluate(forward_point, to_optimize, t_int, LSODA_timeout, Radau_timeout, log10=log10)
            datas.append(forward_data)
            trial_num += 1
            if not success: # increase delta, if evaluation failed
                print(colored(f'\tWarning, forward point failed, {trial_num=}; {current_delta=}', 'yellow'))
                current_delta *= 2
                continue

        # calculate finite difference (gradient[key])
        if success:
            gradient[key] = (forward_data['output'] - central_data['output']) / current_delta
        else:
            gradient[key] = 0.0
    
    return norm_gradient(gradient), datas

def central_difference(point, ranges, to_optimize='energy_efficiency', delta=1e-6, log10=True, t_int=np.array([0.0, 1.0]), LSODA_timeout=30, Radau_timeout=300):
    #keys = ['R_E', 'ratio', 'P_amb', 'alfa_M', 'T_inf', 'surfactant'] + de.excitation_args
    datas = []
    gradient = dict()

    for key in ranges:
        if len(ranges[key]) < 2:
            continue

        # forward point
        forward_point = de.copy(point)
        success_forward = False
        forward_delta = delta
        trial_num = 0
        while not success_forward and trial_num < 3:
            value = point[key] + forward_delta * abs(ranges[key][1] - ranges[key][0])
            forward_point[key] = value
            forward_data, success_forward = evaluate(forward_point, to_optimize, t_int, LSODA_timeout, Radau_timeout, log10=log10)
            trial_num += 1
            datas.append(forward_data)
            if not success_forward: # increase delta, if evaluation failed
                print(colored(f'\tWarning, forward point failed, {trial_num=}; {forward_delta=}', 'yellow'))
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
            backward_data, success_backward = evaluate(backward_point, to_optimize, t_int, LSODA_timeout, Radau_timeout, log10=log10)
            trial_num += 1
            datas.append(backward_data)
            if not success_backward: # increase delta, if evaluation failed
                print(colored(f'\tWarning, backward point failed, {trial_num=}; {backward_delta=}', 'yellow'))
                backward_delta *= 2
                continue

        # central difference
        if success_forward and success_backward:
            gradient[key] = (forward_data['output'] - backward_data['output']) / (forward_delta + backward_delta)
        else:
            print(colored(f'\tError, central difference failed', 'red'))
    
    return norm_gradient(gradient), datas
    
"""________________________________Search________________________________"""

def search(kwargs):
    return gradient_descent(**kwargs)

# TODO describtions
# TODO verbose
# TODO remove log10
# TODO make 2D convergence plot
# TODO remove todos
def gradient_descent(ranges, to_optimize, start_point, step_limit=100, first_step=0.2, decay=0.5, min_step=1e-4, delta=1e-6, verbose=True, log10=True, t_int=np.array([0.0, 1.0]), LSODA_timeout=30, Radau_timeout=300):
    # setup
    solver_kwargs = dict(t_int=t_int, LSODA_timeout=LSODA_timeout, Radau_timeout=Radau_timeout)  # arguments for de.solve()
    keys = [key for key in ranges if not len(ranges[key]) < 2]  # keys of the parameters we want to optimize

    start = time.time()
    step_size = first_step
    step_num = 0
    steps_since_last_decay = 0
    fail_num = 0
    start_point = de.dotdict(start_point)
    
    absolute_best = 1e30
    last_best = 1e30
    datas = []
    last_bests = []

    # learning
    while step_num < step_limit and min_step < step_size and fail_num < 5:
        # calculate gradient (try forward difference, if it fails, try central difference, if that fails too, repeat last step)
        gradient, current_datas = forward_difference(start_point, ranges, to_optimize, delta=delta, **solver_kwargs, log10=log10)
        if gradient is None:
            gradient, current_datas = central_difference(start_point, ranges, to_optimize, delta=delta, **solver_kwargs, log10=log10)
        if gradient is None:
            print(colored(f'\tError, gradient can not be calculated in point {start_point}', 'red'))
            if step_num == 0:
                return [[]], [1.0e30], time.time()-start
            else:
                for key in keys:
                    start_point[key] += change[key]
                start_point = squeeze_into_ranges(start_point, ranges, padding=10.0*delta)
                step_num += 1
                failnum += 1
                continue

        # calculate last best and print stuff
        last_best = min([data['output'] for data in current_datas])
        last_bests.append(last_best)
        absolute_best = min(absolute_best, last_best)
        current_datas = [dict(data) for data in current_datas]
        if verbose:
            print(colored(f'{step_num}. step; {last_best=: .3e}; {absolute_best=: .3e}; {step_size=: .4f}', 'green'))
            point_str = ''.join([f'{key}={start_point[key]: e}; ' for key in keys if isinstance(start_point[key], float)])
            print(f'\tpoint   =({point_str})')

        # calculate change (step's direction)
        change = dict()
        for key in keys:
            interval_width = abs(ranges[key][1] - ranges[key][0])
            change[key] = -step_size * gradient[key] * interval_width

        # decay
        # calculate where the next point would be
        trial_point = de.copy(start_point)
        for key in keys:
            trial_point[key] += change[key]
        trial_point = squeeze_into_ranges(trial_point, ranges, padding=10.0*delta)
        data, success = evaluate(trial_point, to_optimize, **solver_kwargs, log10=log10)
        next_output = data['output']
        current_datas.append(dict(data))
        if success:
            # decay if last decay was too long ago, can't be first decay this way: avoid back and forth steps
            max_step_until_decay = ((1 / decay) + 1 - (1 / decay)%1)
            forced_decay = False
            if step_size != first_step and steps_since_last_decay > max_step_until_decay:
                forced_decay = True
                if verbose: print(colored(f'\tforced decay: {forced_decay}', 'magenta'))
            else:
                forced_decay = False
            # decay if next step is worse than the last best
            if next_output > last_best or forced_decay:
                steps_since_last_decay = 0
                step_size *= decay
                for key in keys:
                    change[key] *= decay
                if verbose: print(colored(f'\tdecayed: {step_size=: .4f}; next_output={next_output}', 'magenta'))
        else:
            if verbose: print(colored(f'\tError, simulation failed in next point', 'red'))

        # make the step
        for key in keys:
            start_point[key] += change[key]
        start_point = squeeze_into_ranges(start_point, ranges, padding=10.0*delta)
        if verbose:
            data, success = evaluate(start_point, to_optimize, **solver_kwargs, log10=log10)
            print(f'\toutput  ={data["output"]}; success={success}')
        
        # print stuff
        datas.append(current_datas)
        step_num += 1
        steps_since_last_decay += 1
        fail_num = 0
        if verbose:
            gradient_str = ''.join([f'{key}={gradient[key]: e}; ' for key in keys])
            print(f'\tgradient=({gradient_str})')
            change_str = ''.join([f'{key}={change[key]: e}; ' for key in keys])
            print(f'\tchange  =({change_str})')

    end = time.time()
    elapsed_time = end-start
    return datas, last_bests, elapsed_time
