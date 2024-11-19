"""Create parameters.py from a .inp and load it"""
# This cell is facultative, you can use an existing parameters.py

# Directory for .inp file:
path = 'INP file examples\\chem_Otomo2018_without_O_reactions_FIXED_by_Cantera.inp'

# import libraries:
import importlib   # for reloading your own files
from termcolor import colored   # for colored error messages

# import inp_data_extractor.py as inp:
try:
    import inp_data_extractor as inp
except ImportError:
    try:
        from Bubble_dynamics_simulation import inp_data_extractor as inp
    except ImportError as error:
        print(colored(f'Error, \'inp_data_extractor.py\' not found', 'red'))
        raise error
    except Exception as error:
        print(colored(f'Error, \'inp_data_extractor.py\' failed to load', 'red'))
        raise error
except Exception as error:
    print(colored(f'Error, \'inp_data_extractor.py\' failed to load', 'red'))
    raise error
importlib.reload(inp)   # reload changes you made

# create parameters.py
inp.extract(path)

# load parameters.py
import parameters as par
importlib.reload(par)
print(par.model)

"""Libraries"""

# for plotting:
#%matplotlib inline
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 18})

import numpy as np   # matrices, math
import time   # runtime measurement
from multiprocessing import Pool, cpu_count   # multithreading
import importlib   # reload changes you made
import itertools   # assemble all combinations of control parameters
import json   # convert dictionary to string

# import full_bubble_model.py as de:
_already_imported = 'de' in globals()
try:
    import full_bubble_model as de
except ImportError:
    try:
        from  Bubble_dynamics_simulation import full_bubble_model as de
    except ImportError as _error:
        print(colored(f'Error, \'full_bubble_model.py\' not found', 'red'))
        raise _error
    except Exception as _error:
        print(colored(f'Error, \'full_bubble_model.py\' failed to load', 'red'))
        raise _error
except Exception as _error:
    print(colored(f'Error, \'full_bubble_model.py\' failed to load', 'red'))
    raise _error
if _already_imported: importlib.reload(de)   # reload changes you made

"""Control parameter ranges and division"""
# a list for each control parameter, containing all the possible values

ranges = dict(
  # Initial conditions:
    # bubble equilibrium radius [m]
    R_E = [1.0e-6*x for x in [100.0]],#[1.0e-6*x for x in np.logspace(0.0,4.0,num=201)],#
    # initial radius / equilibrium radius R_0/R_E [-]
    ratio = [x for x in [1.0]],#[x for x in np.linspace(1.2,20,num=377)],#
  # Ambient parameters:
    # Standard pressure [Pa]
    P_standard=[101325.0], # [Pa]
    # ambient pressure [Pa]
    P_amb = [x * par.bar2Pa for x in np.logspace(0.0,3.0,2)], # [bar --> Pa]
    # ambient temperature [K]       
    T_inf = [par.absolute_zero + x for x in [20.0]], #par.absolute_zero +  # [°C --> K]
    # Molar fractions of species in the initial bubble (H2 and N2) [-]
    fractions = [[x,1.0-x] for x in [0.75]],# in np.linspace(0.4,0.85,num=10)],
  # Liquid parameters:
    # Surface tension modifier [-]
    surfactant = [1.0],
    # water accommodation coefficient [-]
    alfa_M = [0.35], #[-],
    Gamma=[1.0], #[-],
    sigma_evap=[0.4], #[-],
    C_4_starred=[-2.1412], #[-]
    #Dynamic viscosity of the liquid [Pa*s]
    mu_L = [0.001 * x for x in [1.0]],#,
    #Sound velocity in the liquid [m/s]
    c_L = [1483.0 * x for x  in [1.0]],#,
    # vapour pressure [Pa]
    #P_v = #par.P_v, # calculated from T_inf
    # dynamic viscosity [Pa*s]
    #mu_L = [par.mu_L], # calculated from T_inf
    # density [kg/m^3]
    rho_L =  [998.20],
    # sound speed [m/s]
    #c_L = [par.c_L],
  # Excitation parameters: (excitation_type = sin_impulse)
    # excitation amplitude [Pa]
    p_A = [-x * par.bar2Pa for x in [2.0]],#np.linspace(1.0, 3.0, 50)], # [bar --> Pa]
    # excitation frequency [Hz]
    freq =  [20000.00],
    #2-frequency case:
    #Exciting frequency 1 [Hz]
    freq1=[2.0e4],                          # [Hz]
    #Frequency ratio:
    freq_ratio=[2.0],#freq2=1.0e5,            freq2/freq1  # [-]
    #Exciting pressure amplitude 1 [Pa]
    p_A1=[-2.0e5],                            # [Pa]
    #Exciting pressure amplitude 2 [Pa]
    p_A2=[-0.8e5],                            # [Pa]
    # excitation duration in period times [-]
    n =  [1.00],
    
    #Thermodynamical case: Constant volume...
    thermodynamicalcase = [2], #0 = 'NonIsothermal-ConstantVolume'
)

for key in de.excitation_args:
    if key not in ranges:
        print(colored(f'Error, {key} not in ranges', 'red'))
# print total combinations:
for key in ranges:
    print(f'{key}: {len(ranges[key])}')
total_combinations = f'total combinations: {np.prod([len(ranges[key]) for key in ranges])}'
print(''.join(['_' for i in range(len(total_combinations))]))
print(total_combinations)

"""Get all combinations"""
# Make a list, with one dictionary for eachy parameter combinations

start = time.time()
cpars = []
ID = 1
for values in itertools.product(*ranges.values()):
    cpar = dict(zip(ranges.keys(), values))
    cpar['ID'] = ID                      # ID of control parameter (not used during calculation)
    cpar['gases'] = [par.index['H2'], par.index['N2']]    # indexes of species in initial bubble (list of species indexes)
    cpar['fractions'] = [0.75, 0.25]            # molar fractions of species in initial bubble (list of fractions for every gas)
    # Calculate pressure/temperature dependent parameters:
    cpar['mu_L'] = de.viscosity(cpar['T_inf'])
    cpar['P_v'] = de.vapour_pressure(cpar['T_inf'])
    cpar['theta_phase']=3.0*np.pi/2.0*(1.0-cpar['freq_ratio'])
    cpars.append(cpar)
    ID += 1

print(f'Assemble cpars: {time.time()-start:.2f} s')
start = time.time()

# Create input dictionary for de.simulate(), a list of dictionaries with cpar and other arguments
kwargs_list = [dict(cpar=cpar, t_int=np.array([0.0, 1.0e4]), LSODA_timeout=30, Radau_timeout=100) for cpar in cpars]
end = time.time()
print(f'Assemble kwargs_list: {time.time()-start:.2f} s')

"""Save settings as txt"""

# create folder for parameter study results:
file = de.Make_dir('test_1atm_20000Hz_2D')

# save all settings (full_bubble_model.py, parameters.py, ranges) as txt:
combined_str = f'''
total_combinations = {np.prod([len(ranges[key]) for key in ranges])}
ranges = {json.dumps(ranges, indent=4)}'''

file.write_string(combined_str, 'brutefroce_parameter_sweep_settings')

"""Parameter study, multithread"""
# Runs each combinations in cpars, and saves the results into CSV files
# use Pool(processes=cpu_count()-1) to limit number of threads being used.
# use pool.imap(...) instead of pool.imap_unordered(...) to preserve order in which cpars was made

max_lines = 10000    # maximum length of a CSV
best_energy_demand = 1e30

start = time.time()
file.new_file()
with Pool(processes=cpu_count(), maxtasksperchild=1000) as pool:
    results = pool.imap_unordered(de.simulate, kwargs_list)

    for data in results:
      # save results:
        if file.lines > max_lines:
            file.close()
            file.new_file()
        data = de.dotdict(data)
        file.write_line(data)
      # print stuff:
        if data.energy_demand > 0 and data.energy_demand < best_energy_demand:
            best_energy_demand = data.energy_demand
        excitation_params = ''.join([f'{key}={data[key]: <12}; ' for key in de.excitation_args])
        print(f'index: {data.ID: >8}/{len(cpars)};   error_code: {data.error_code: >4} (success={data.success});   steps: {data.steps: <8};   runtime: {data.elapsed_time: 6.2f} [s]   |   '+
              f'R_E={1e6*data.R_E: 6.2f} [um]; ratio={data.ratio: 6.2f} [-]; P_amb={1e-5*data.P_amb: 6.2f} [bar]; T_inf={data.T_inf-273.15: 6.2f} [°C]; '+
              f'alfa_M={data.alfa_M: 6.2f} [-]; P_v={data.P_v: 6.0f} [Pa]; mu_L={data.mu_L: 6.5f} [Pa*s]; rho_L={data.rho_L: 6.1f} [kg/m^3]; c_L={data.c_L: 6.0f} [m/s]; '+
              f'surfactant={100*data.surfactant: 3.0f} [%]   |   {excitation_params}   |   '+
              f'{de.target_specie} production: {data.energy_demand: e} [MJ/kg] (best: {best_energy_demand: .1f} [MJ/kg])'+
              '                                                 ', end='\r')
              
file.close()
end = time.time()
elapsed = end - start
print(f'\n\nDONE')
print(f'total time: {(elapsed / 3600): .0f} hours {((elapsed % 3600) / 60): .0f} mins')
print(f'            {elapsed: .2f} [s]   ({(elapsed / len(    cpars)): .2f} [s/run])')