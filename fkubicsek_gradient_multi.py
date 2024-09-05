"""Create parameters.py from a .inp and load it"""
# This cell is facultative, you can use an existing parameters.py

# Directory for .inp file:
path = 'INP file examples//chem_Otomo2018_without_O_reactions_FIXED_by_Cantera.inp'#'INP file examples//chem_Otomo2018.inp'#'INP file examples//chemkin_AR_HE.inp'

# import libraries:
import math
import importlib   # for reloading your own files
import itertools   # assemble all combinations of control parameters
import json   # convert dictionary to string

# my own files:
import parameters as par

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


# my own files:
try:
    import inp_data_extractor as inp
except:
    try:
        import Bubble_dynamics_simulation.inp_data_extractor as inp
    except:
        print(colored(f'Error, \'inp_data_extractor.py\' not found', 'red'))
importlib.reload(inp)

# create parameters.py
inp.extract(path)

# load parameters.py
import parameters as par
importlib.reload(par)
print(par.model)

"""Libraries"""

# for plotting:
#%matplotlib inline 
#If you run this code on the HDS department CPU you should comment the previous line! This line is needed only if you run it in Jupyter IPython.
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 18})

import numpy as np   # matrices, math
import time   # runtime measurement
import random   # random number generator
from multiprocessing import Pool, cpu_count   # multithreading
import importlib   # reload changes you made

# my own file:
import gradient_method as gm    # full bubble model
importlib.reload(gm)   # reload changes you made

to_optimize = 'energy_efficiency'   # key in data from de.get_data()

"""Starting points and control parameter ranges"""
# Starting points
starting_points = dict(
    ID = [0], #default value, will be overwritten
    # Equilibrium radius [um --> m]
    R_E = [1.0e-6*x for x in [1.0e6*0.00966115244195245]], #[1.0, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0, 15.0, 20.0, 25.0]],
    # R_star / R_E [-]
    ratio = [5.74579143637999], #[2.0, 3.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0],
    # Ambient pressure [Pa]
    P_amb = [par.bar2Pa*x for x in [1000.0]],#[0.00, 0.50]], #[1e5*x for x in [1.0, 10.0, 25.0, 50.0]],
    # Accommodation coeff. for evaporation [-]
    alfa_M = [0.35],
    # Ambient temperature [°C --> K]
    T_inf = [273.15+x for x in [20.0]],
    # Surface tension modifier [-]
    surfactant = [1.0],#[0.0,1.0],
    # Molar fractions of species in the initial bubble (H2 and N2) [-]
    fractions = [[x,1.0-x] for x in [0.7302360446128061]], 
    #Exciting frequency 1 [Hz]
    freq1 = [2.0e4], 
    #Exciting frequency 2 [Hz]
    freq2 = [1.0e5],
    #Exciting pressure amplitude 1 [Pa]
    pA1 = [0.0e5], 
    #Exciting pressure amplitude 2 [Pa]
    pA2 = [0.0e5], 
    #Initial phase angle between the exciting pressure waves [rad]
    theta_phase = [0.1],
    # excitation amplitude [Pa]
    p_A = [-x * par.bar2Pa for x in [0.0]],#np.linspace(1.0, 3.0, 50)], # [bar --> Pa]
    # excitation frequency [Hz]
    freq =  [20000.00],
    # excitation duration in period times [-]
    n =  [1.00],
    #Dynamic viscosity of the liquid [Pa*s]
    mu_L = [0.001 * x for x in np.logspace(np.log10(0.2),np.log10(100.0),num=8)],#[0.1,100.0]],
    #Sound velocity in the liquid [m/s]
    c_L = [1483.0 * x for x  in np.linspace(0.7,1.3,num=8)],#[0.0,2.0]], 
    #Density [kg/m^3]
    rho_L =  [998.20],
    #Thermodynamical case: Constant volume...
    thermodynamicalcase = [2], #0 = 'NonIsothermal-ConstantVolume', 1 = 'Isothermal-ConstantVolume', other: Keller-Miksis
)

ID=1
# Range or a single value for most control parameters
ranges = dict(
    # Equilibrium radius [um --> m]
    R_E = [1.0e-6*x for x in [9000.0,10000.0]], #[1.0, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0, 15.0, 20.0, 25.0]],
    # R_star / R_E [-]
    ratio = [1.1,20.0], #[2.0, 3.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0],
    # Ambient pressure [Pa]
    P_amb = [par.bar2Pa*x for x in [1000.0]],#[0.00, 0.50]], #[1e5*x for x in [1.0, 10.0, 25.0, 50.0]],
    # Accommodation coeff. for evaporation [-]
    alfa_M = [0.35],
    # Ambient temperature [°C --> K]
    T_inf = [273.15+x for x in [20.0]],
    # Surface tension modifier [-]
    surfactant = [1.0],#[0.0,1.0],
    # Molar fractions of species in the initial bubble (H2 and N2) [-]
    fractions = [[0.40, 0.60],[0.85, 0.15]], 
    #Exciting frequency 1 [Hz]
    freq1 = [2.0e4], 
    #Exciting frequency 2 [Hz]
    freq2 = [1.0e5], 
    #Exciting pressure amplitude 1 [Pa]
    pA1 = [0.0e5], 
    #Exciting pressure amplitude 2 [Pa]
    pA2 = [0.0e5], 
    #Initial phase angle between the exciting pressure waves [rad]
    theta_phase = [0.1],
    # excitation amplitude [Pa]
    p_A = [-x * par.bar2Pa for x in [0.0]],#np.linspace(1.0, 3.0, 50)], # [bar --> Pa]
    # excitation frequency [Hz]
    freq =  [20000.00],
    # excitation duration in period times [-]
    n =  [1.00],
    #Dynamic viscosity of the liquid [Pa*s]
    mu_L = [0.001 * x for x in [1.0]],#[0.1,100.0]],
    #Sound velocity in the liquid [m/s]
    c_L = [1483.0 * x for x  in [1.0]],#[0.0,2.0]], 
    #Density [kg/m^3]
    rho_L =  [998.20],
)

""" Gradient method, multithread"""
searches = 1    # number os total searches
for key in starting_points:
    searches *= len(starting_points[key])

best_output = 1.0e30
total_point_num = 0
num = 0
to_plot = []
last_points = []
start = time.time()
file = gm.de.Make_dir('test 1')

cpars = []
ID = 0
for values in itertools.product(*starting_points.values()):
    cpar = dict(zip(starting_points.keys(), values))
    cpar['ID'] = ID                      # ID of control parameter (not used during calculation)
    cpar['gases'] = [par.index['H2'], par.index['N2']]    # indexes of species in initial bubble (list of species indexes)
    if not 'fractions' in cpar:
        cpar['fractions'] = [0.75, 0.25]            # molar fractions of species in initial bubble (list of fractions for every gas)
    # Calculate pressure/temperature dependent parameters:
    if not 'mu_L' in cpar:
        cpar['mu_L'] = de.viscosity(cpar['T_inf'])
    if not 'P_v' in cpar:
        cpar['P_v'] = de.vapour_pressure(cpar['T_inf'])
    cpars.append(cpar)
    ID += 1

kwargs_list = [dict(
    ranges=ranges,
    to_optimize=to_optimize,
    start_point=gm.fix_starter_point(ranges, i, cpars[i], gases=[gm.de.par.index['H2'], gm.de.par.index['N2']]),
    max_steps=2000, #in gradient method
    first_step=0.05, #between two parameter combinations
    min_step=0.00001, #between two parameter combinations
    decay=0.6,
    gamma=1.0, #x_n+1=x_n-gamma*grad(F(x_n))
    t_int=np.array([0.0, 1.0e4]),
    LSODA_timeout=30,
    Radau_timeout=300,
    ) for i in range(ID)]

# save points
file.new_file()
with Pool(processes=2, maxtasksperchild=10) as pool: #processes=cpu_count()
    results = pool.imap_unordered(gm.search,kwargs_list)
    for result in results:
        datas, best_central_values, elapsed, step_num, point_num = result
        total_point_num += point_num
        num += 1
        to_plot.append(best_central_values)
        last_points.append(best_central_values[-1])
        for i in range(len(datas)):
            if datas[i]['output'] < best_output:
                best_output = datas[i]['output']

        for data in datas:
            file.write_line(data)
        
        # print stuff:
        print(f'{num}/{searches}: Total {step_num} steps and {point_num} points, finished in {elapsed: .2f} [s]   ({(elapsed / point_num): .2f} [s/run]).   '+
              f'Final {to_optimize} (central value): {best_central_values[-1]: .1f} (best: {best_output: .1f})')
            
file.close()
end = time.time()
elapsed = end - start
print(f'\n\nDONE')
print(f'total time: {((elapsed-elapsed % 3600) / 3600): .0f} hours {((elapsed % 3600) / 60): .0f} mins')
print(f'            {elapsed: .2f} [s]   ({(elapsed / total_point_num): .2f} [s/run])')

"""Plot convergence of last searches"""

plt.rcParams.update({'font.size': 16})
if len(to_plot) > 10:
    to_plot2 = to_plot[-10:]
else:
    to_plot2 = to_plot

fig, ax = plt.subplots(1, 1, figsize=(15, 6))
fig.suptitle('Convergence of last 10 searches', fontsize=22)
ax.set_ylabel(f'{to_optimize}')
ax.set_ylim(100.0, 0.5e5)
ax.set_yscale('log')
ax.set_xlabel('steps [-]')
ax.grid()
for plot in to_plot2:
    ax.plot(plot, '.-', linewidth=1.0)
plt.show()

"""Plot the distribution of the optimums"""

plt.rcParams.update({'font.size': 18})
ranges2 = dict()
for key in ranges:
    if len(ranges[key]) > 1:
        ranges2[key] = ranges[key]

n = len(ranges2)
fig, ax = plt.subplots(n, 1, figsize=(24, 5*n))
for i, key in enumerate(ranges2):
    last_values = [last_point[key] for last_point in last_points]
    y = [i for i, last_point in enumerate(last_points)]
        
  # ploting:
    ax[i].scatter(last_values, y, s=25, color='b')
    ax[i].set_ylabel('ID of search')
    ax[i].set_xlabel(key)
    ax[i].set_xlim(ranges2[key])
    n = len(last_points)
    ax[i].set_ylim([-0.5*n, 1.5*n])
    ax[i].grid()
    # best point:
    index = [x[-1] for x in to_plot].index(best_output)
    best_point = last_points[index]
    ax[i].scatter([best_point[key]], [0.5*n], s=200, color='r')