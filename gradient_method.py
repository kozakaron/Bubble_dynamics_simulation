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
def rand_point(ranges, ID=0, gases=[de.par.index['AR']], first_step=0.1):
    point = dict(ID=ID)
    for key in ranges:
        if len(ranges[key]) < 2:
            point[key] = ranges[key][0]
        elif (key == 'fractions'):
            x=random.uniform(ranges[key][0][0] + first_step*abs(ranges[key][1][0] - ranges[key][0][0]), ranges[key][1][0] - first_step*abs(ranges[key][1][0] - ranges[key][0][0])) #fraction of the first material (in case of NH3 it is H2)
            point[key] = [x,1.0-x]
        else:
            point[key] = random.uniform(ranges[key][0] + first_step*abs(ranges[key][1] - ranges[key][0]), ranges[key][1] - first_step*abs(ranges[key][1] - ranges[key][0]))
    point['gases'] = gases
    point['P_v'] = de.VapourPressure(point['T_inf'])  # [Pa]
    if not 'mu_L' in point:
        point['mu_L'] = de.Viscosity(point['T_inf'])  # [Pa*s]
    print(point)
    return point

# generates a cpar dict (point) with specified values
def fix_starter_point(ranges, ID=0, gases=[de.par.index['AR']], first_step=0.1):
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
    cpar = de.dotdict(point)
    cpar.P_v = de.VapourPressure(T=cpar.T_inf) # [Pa]
    if not 'mu_L' in cpar:
        cpar.mu_L = de.Viscosity(T=cpar.T_inf) # [Pa]
    num_sol, error_code, elapsed_time, negativepressure_error = de.solve(cpar, t_int, LSODA_timeout, Radau_timeout)   # simulation without plotting
    data = de.get_data(cpar, num_sol, error_code, elapsed_time, negativepressure_error)   # post processing
    # return value to optimize
    output = data[to_optimize]
    if output < 0.0 and data['expansion_work'] > 0.0:
        output = 1.0e30
    if data['error_code'] > 3 or negativepressure_error:
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
    return gradient_method(**kwargs)

def gradient_method(ranges, to_optimize, start_point, max_steps=100, first_step=0.1, min_step=0.02, decay=0.6, gamma=1.0, t_int=np.array([0.0, 1.0]), LSODA_timeout=30, Radau_timeout=300):
    start = time.time()
    point_num = 0
    step_num = 0
    current_step = first_step
    start_point['step_size'] = current_step
    data = evaluate(start_point, to_optimize, t_int, LSODA_timeout, Radau_timeout)
    best_central_points = [start_point]
    datas = [copy_dotdict(data)]
    keys = ['R_E', 'ratio', 'P_amb', 'alfa_M', 'T_inf', 'surfactant', 'fractions', 'freq1', 'freq2', 'pA1', 'pA2', 'theta_phase', 'mu_L', 'c_L']
    
    counter = 0
    while step_num < max_steps and current_step > min_step:
        best_central_points[-1]['step_size'] = current_step
        new_points = []   
        difference = np.zeros(len(keys), dtype=np.float64)
        step_size = np.zeros(len(keys), dtype=np.float64)
        notinrange = 0
        
        for i in range(len(keys)):   # make two steps in all directions 
            key = keys[i]
            if len(ranges[key]) < 2:
                continue
            if key != 'fractions':
                step_size[i] = current_step * abs(ranges[key][1] - ranges[key][0])

                #step 1
                new_point = copy_dotdict(best_central_points[-1])
                new_point[key] -= step_size[i]
                if in_range(new_point, ranges):
                    data = evaluate(new_point, to_optimize, t_int, LSODA_timeout, Radau_timeout)
                    new_points.append(new_point)
                    datas.append(copy_dotdict(data))
                    point_num += 1
                else:
                    print('Warning! Step_num: '+str(step_num)+', the new point is not in range!')
                    print('key:'+str(key));print('ranges:');print(ranges);print('new_point:');print(new_point)
                    notinrange = 1

                #step 2
                new_point = copy_dotdict(best_central_points[-1])
                new_point[key] += step_size[i]
                if in_range(new_point, ranges):
                    data = evaluate(new_point, to_optimize, t_int, LSODA_timeout, Radau_timeout)
                    new_points.append(new_point)
                    datas.append(copy_dotdict(data))
                    point_num += 1
                    difference[i] = (datas[-1]['output']-datas[-2]['output']) / 2.0
                else:
                    print('Warning! Step_num: '+str(step_num)+', the new point is not in range!')
                    print('key:'+str(key));print('ranges:');print(ranges);print('new_point:');print(new_point)
                    notinrange = 1
                
            else:
                step_size[i] = current_step * abs(ranges[key][1][0] - ranges[key][0][0])

                #step 1
                new_point = copy_dotdict(best_central_points[-1])
                new_point[key][0] -= step_size[i]
                new_point[key][1] += step_size[i]
                if in_range(new_point, ranges):
                    data = evaluate(new_point, to_optimize, t_int, LSODA_timeout, Radau_timeout)
                    new_points.append(new_point)
                    datas.append(copy_dotdict(data))
                    point_num += 1
                else:
                    print('Warning! Step_num: '+str(step_num)+', the new point is not in range!')
                    print('key:'+str(key));print('ranges:');print(ranges);print('new_point:');print(new_point)
                    notinrange = 1
                    difference[i] = 0.0

                #step 2
                new_point = copy_dotdict(best_central_points[-1])
                new_point[key][0] += step_size[i]
                new_point[key][1] -= step_size[i]
                if in_range(new_point, ranges):
                    data = evaluate(new_point, to_optimize, t_int, LSODA_timeout, Radau_timeout)
                    new_points.append(new_point)
                    datas.append(copy_dotdict(data))
                    point_num += 1
                    difference[i] = (datas[-1]['output']-datas[-2]['output']) / 2.0
                else:
                    print('Warning! Step_num: '+str(step_num)+', the new point is not in range!')
                    print('key:'+str(key));print('ranges:');print(ranges);print('new_point:');print(new_point)
                    notinrange = 1
                    difference[i] = 0.0
            
        step_num += 1;print('step_num: '+str(step_num)+', max_steps: '+str(max_steps))

        # determine best step
        if step_num == 0:
            continue
        values = [new_point['output'] for new_point in new_points]
        
        if notinrange == 0:
            best_value = min(values)
            best_point = new_points[values.index(best_value)]
            
            # if it's better than the previous then next step
            if best_value < best_central_points[-1]['output']:
                length_difference = 0.0
                for i in range(len(keys)):
                    length_difference += difference[i]**2.0
                length_difference = length_difference**0.5
                if length_difference == 0.0:
                    length_difference = 1.0; current_step *= decay
                re_length_difference = 1.0 / length_difference

                for i in range(len(keys)):
                    key = keys[i]
                    if key != 'fractions':
                        best_point[key] = (best_central_points[-1][key] - gamma*difference[i]*re_length_difference*step_size[i]) 
                    else:
                        best_point[key][0] = (best_central_points[-1][key][0] - gamma*difference[i]*re_length_difference*step_size[i])
                        best_point[key][1] = (best_central_points[-1][key][1] + gamma*difference[i]*re_length_difference*step_size[i]) 
                best_central_points.append(best_point)
                data = evaluate(best_point, to_optimize, t_int, LSODA_timeout, Radau_timeout)
                datas.append(copy_dotdict(data))
                print('best central point: ');print(best_point)
                counter += 1
                if counter > 19:
                    current_step *= decay
                    counter = 0
            else:   # if new step is worse then retake the step with smaller step_size
                current_step *= decay
        else:   # if we are not in range then retake the step with smaller step_size
            current_step *= decay

    end = time.time()
    elapsed = end - start
    return [datas, [point['output'] for point in best_central_points], elapsed, step_num, point_num]
