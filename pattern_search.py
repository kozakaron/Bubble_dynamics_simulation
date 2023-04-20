"""________________________________Libraries________________________________"""

import numpy as np   # matrices, math
import time   # runtime measurement
import random   # random number generator
import importlib   # reload changes you made

# my own file:
import full_bubble_model as de    # full bubble model
importlib.reload(de)   # reload changes you made

"""________________________________Helper functions________________________________"""

# generates a random cpar dict (point) within ranges
def rand_point(ranges, ID=0, gases=[de.par.index['AR']], fractions=[1.0]):
    point = dict(ID=ID)
    for key in ranges:
        if len(ranges[key]) < 2:
            point[key] = ranges[key][0]
        else:
            point[key] = random.uniform(ranges[key][0], ranges[key][1])
    point['gases'] = gases
    point['fractions'] = fractions
    point['P_v'] = de.VapourPressure(point['T_inf'])  # [Pa]
    point['mu_L'] = de.Viscosity(point['T_inf'])  # [Pa*s]
    return point

# checks if a cpar is in ranges or not
def in_range(point, ranges):
    for key in ranges:
        if key in ['ID', 'gases', 'fractions']:
            continue
        if len(ranges[key]) < 2:
            if point[key] != ranges[key][0]:
                return False
        elif point[key] < ranges[key][0] or ranges[key][1] < point[key]:
            return False
    return True

# runs the simulation, and returns with data extended with output (what we want to optimize)
def evalvuate(point, to_optimize, t_int, LSODA_timeout, Radau_timeout):
    cpar = de.dotdict(point)
    cpar.P_v = de.VapourPressure(T=cpar.T_inf) # [Pa]
    cpar.mu_L = de.Viscosity(T=cpar.T_inf) # [Pa]
    num_sol, error_code, elapsed_time = de.solve(cpar, t_int, LSODA_timeout, Radau_timeout)   # simulation without plotting
    data = de.get_data(point, num_sol, error_code, elapsed_time)   # post processing
    
    # return value to optimize
    output = data[to_optimize]
    if output < 0.0:
        output = 1.0e30
    if data['error_code'] > 3:
        output = 1.0e30
    point['output'] = output
    data['output'] = output
    return data

# makes a deep copy of a dotdict
# normally python only copies pointers
def copy_dotdict(inp):
    ret = dict()
    for key in inp:
        if type(inp[key]) == list or type(inp[key]) == np.ndarray:
            ret[key] = [element for element in inp[key]]
        else:
            ret[key] = inp[key]
    return ret


"""________________________________Search________________________________"""

def search(kwargs):
    return pattern_search(**kwargs)

def pattern_search(ranges, to_optimize, start_point, max_steps=100, first_step=0.1, min_step=0.02, decay=0.6, t_int=np.array([0.0, 1.0]), LSODA_timeout=30, Radau_timeout=300):
    start = time.time()
    point_num = 0
    step_num = 0
    current_step = first_step
    start_point['step_size'] = current_step
    data = evalvuate(start_point, to_optimize, t_int, LSODA_timeout, Radau_timeout)
    best_points = [start_point]
    datas = [copy_dotdict(data)]

    while step_num < max_steps and current_step > min_step:
        best_points[-1]['step_size'] = current_step
        new_points = []
        for key in ['R_E', 'ratio', 'P_inf', 'alfa_M', 'T_inf', 'surfactant']:   # make 2 steps in all directions
            if len(ranges[key]) < 2:
                continue
            step_size = current_step * abs(ranges[key][1] - ranges[key][0])

            # step 1
            new_point = copy_dotdict(best_points[-1])
            new_point[key] -= step_size
            if in_range(new_point, ranges):
                data = evalvuate(new_point, to_optimize, t_int, LSODA_timeout, Radau_timeout)
                new_points.append(new_point)
                datas.append(copy_dotdict(data))
                point_num += 1

            # step 2
            new_point = copy_dotdict(best_points[-1])
            new_point[key] += step_size
            if in_range(new_point, ranges):
                data = evalvuate(new_point, to_optimize, t_int, LSODA_timeout, Radau_timeout)
                new_points.append(new_point)
                datas.append(copy_dotdict(data))
                point_num += 1

        step_num += 1

        # determine best step
        if step_num == 0:
            continue
        values = [new_point['output'] for new_point in new_points]
        best_value = min(values)
        best_point = new_points[values.index(best_value)]
        # if it's better than the previous than next step
        if best_value < best_points[-1]['output']:
            best_points.append(best_point)
        else:   # if new step is worse than retake the step with smaller step_size
            current_step *= decay

    end = time.time()
    elapsed = end - start
    return [datas, [point['output'] for point in best_points], elapsed, step_num, point_num]