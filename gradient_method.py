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

# generates a random cpar dict (point) within ranges
def rand_point(ranges, gases, ID=0, padding=0.1):
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

    point['gases'] = gases
    if not 'P_v' in point:
        point['P_v'] = de.VapourPressure(point['T_inf'])  # [Pa]
    if not 'mu_L' in point:
        point['mu_L'] = de.Viscosity(point['T_inf'])  # [Pa*s]

    return point

# generates a cpar dict (point) with specified values
def fix_starter_point(ranges, ID=0, gases=[de.par.index['H2']], first_step=0.1):
    point = dict(ID=ID)
    for key in ranges:
        if len(ranges[key]) < 2:
            point[key] = ranges[key][0]
    point['R_E'] = 40.0e-6
    point['ratio'] = 6.2
    point['P_amb'] = 111457.5
    point['fractions'] = [0.70,0.30]
    point['gases'] = gases
    point['P_v'] = de.VapourPressure(point['T_inf'])  # [Pa]
    if not 'mu_L' in point:
        point['mu_L'] = de.Viscosity(point['T_inf'])  # [Pa*s]
    return point

# checks if a cpar is in ranges or not
def in_range(point, ranges):
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

# runs the simulation, and returns with data extended with output (what we want to optimize)
def evaluate(point, to_optimize, t_int, LSODA_timeout, Radau_timeout):
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
    if not success:
        output = 1.0e30
    elif output < 0.0 and data['expansion_work'] > 0.0:
        output = 1.0e30
    #output = np.log10(output)
    point['output'] = output
    data['output'] = output
    ret = {key: de.copy(data[key]) for key in data}
    return ret, success


"""________________________________Search________________________________"""

def search(kwargs):
    return gradient_method(**kwargs)

def gradient_method(ranges, to_optimize, start_point, max_steps=100, first_step=0.1, min_step=0.0002, decay=0.6, learning_rate=1.0, t_int=np.array([0.0, 1.0]), LSODA_timeout=30, Radau_timeout=300):
    start = time.time()
    point_num = 0
    step_num = 0
    current_step = first_step
    start_point['step_size'] = current_step
    best_central_points = [start_point]
    data, success = evaluate(start_point, to_optimize, t_int, LSODA_timeout, Radau_timeout)
    datas = [de.copy(data)]
    keys = ['R_E', 'ratio', 'P_amb', 'alfa_M', 'T_inf', 'surfactant', 'fractions', 'mu_L', 'c_L'] + de.excitation_args
    
    counter = 0
    while step_num < max_steps and current_step > min_step:
        best_central_points[-1]['step_size'] = current_step
        new_points = []
        difference = np.zeros(len(keys), dtype=np.float64)
        step_size = np.zeros(len(keys), dtype=np.float64)
        notinrange = False
        
        for i, key in enumerate(keys):   # make two steps in all directions 
            if len(ranges[key]) < 2:
                continue
            if key != 'fractions':
                step_size[i] = current_step * abs(ranges[key][1] - ranges[key][0])

                #step 1
                new_point = de.copy(best_central_points[-1])
                new_point[key] -= step_size[i]
                if in_range(new_point, ranges):
                    data, success = evaluate(new_point, to_optimize, t_int, LSODA_timeout, Radau_timeout)
                    new_points.append(new_point)
                    datas.append(de.copy(data))
                    point_num += 1
                else:
                    print(colored(f'Warning, in step {step_num=}, the new point is not in range!', 'yellow'))
                    print(colored(f'\t{key=}; ranges={ranges[key]}; current_value={new_point[key]}; {new_point=}', 'yellow'))
                    notinrange = True

                #step 2
                new_point = de.copy(best_central_points[-1])
                new_point[key] += step_size[i]
                if in_range(new_point, ranges):
                    data, success = evaluate(new_point, to_optimize, t_int, LSODA_timeout, Radau_timeout)
                    new_points.append(new_point)
                    datas.append(de.copy(data))
                    point_num += 1
                    difference[i] = (datas[-1]['output'] - datas[-2]['output']) / 2.0
                else:
                    print(colored(f'Warning, in step {step_num=}, the new point is not in range!', 'yellow'))
                    print(colored(f'\t{key=}; ranges={ranges[key]}; current_value={new_point[key]}; {new_point=}', 'yellow'))
                    notinrange = True
                    difference[i] = 0.0
                
            else:
                step_size[i] = current_step * abs(ranges[key][1][0] - ranges[key][0][0])

                #step 1
                new_point = de.copy(best_central_points[-1])
                new_point[key][0] -= step_size[i]
                new_point[key][1] += step_size[i]
                if in_range(new_point, ranges):
                    data, success = evaluate(new_point, to_optimize, t_int, LSODA_timeout, Radau_timeout)
                    new_points.append(new_point)
                    datas.append(de.copy(data))
                    point_num += 1
                else:
                    print(colored(f'Warning, in step {step_num=}, the new point is not in range!', 'yellow'))
                    print(colored(f'\t{key=}; ranges={ranges[key]}; current_value={new_point[key]}; {new_point=}', 'yellow'))
                    notinrange = True

                #step 2
                new_point = de.copy(best_central_points[-1])
                new_point[key][0] += step_size[i]
                new_point[key][1] -= step_size[i]
                if in_range(new_point, ranges):
                    data, success = evaluate(new_point, to_optimize, t_int, LSODA_timeout, Radau_timeout)
                    new_points.append(new_point)
                    datas.append(de.copy(data))
                    point_num += 1
                    difference[i] = (datas[-1]['output'] - datas[-2]['output']) / 2.0
                else:
                    print(colored(f'Warning, in step {step_num=}, the new point is not in range!', 'yellow'))
                    print(colored(f'\t{key=}; ranges={ranges[key]}; current_value={new_point[key]}; {new_point=}', 'yellow'))
                    notinrange = True
                    difference[i] = 0.0
            
        step_num += 1
        print(colored(f'{step_num=}/{max_steps}', 'green'))

        # determine best step
        if step_num == 0:
            continue
        values = [new_point['output'] for new_point in new_points]
        
        if notinrange == False:
            best_value = min(values)
            best_point = new_points[values.index(best_value)]
            
            # if it's better than the previous then next step
            if best_value < best_central_points[-1]['output']:
                length_difference = np.linalg.norm(difference) # length of the difference vector
                if length_difference == 0.0:
                    length_difference = 1.0
                    current_step *= decay
                re_length_difference = 1.0 / length_difference

                for i, key in enumerate(keys):
                    if key != 'fractions':
                        best_point[key] = (best_central_points[-1][key] - learning_rate * difference[i] * re_length_difference * step_size[i]) 
                    else:
                        best_point[key][0] = (best_central_points[-1][key][0] - learning_rate * difference[i] * re_length_difference * step_size[i])
                        best_point[key][1] = (best_central_points[-1][key][1] + learning_rate * difference[i] * re_length_difference * step_size[i])

                # calculate new point
                best_central_points.append(best_point)
                data, success = evaluate(best_point, to_optimize, t_int, LSODA_timeout, Radau_timeout)
                datas.append(de.copy(data))
                print(f'best central point: output={best_point["output"]}; {best_point=}')
                counter += 1
                if counter > 19:
                    current_step *= decay
                    counter = 0
            else:   # if new step is worse then retake the step with smaller step_size
                current_step *= decay
        else:   # if we are not in range then retake the step with smaller step_size
            current_step *= decay
        print(f'{current_step=}; {step_num=}')

    end = time.time()
    elapsed = end - start
    return [datas, [point['output'] for point in best_central_points], elapsed, step_num, point_num]
