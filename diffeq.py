"""________________________________Libraries________________________________"""

# for plotting:
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 18})

import numpy as np   # matrices, math
from scipy.integrate import solve_ivp   # differential equation solver
from scipy.signal import argrelmin   # loc min finding
import time   # runtime measurement
from numba import jit, njit   # Just In Time compiler
from numba.types import unicode_type, float64, float32, int64, int32   # JIT types
import os    # file management
import importlib   # For reloading your own files

# my own files:
import chemkin_AR_HE as par   # numeric constants and coefficents
importlib.reload(par)   # reload changes you made

# dot.notation access to dictionary attributes
# instead of dictionary['key'] you can use dictionary.key
class dotdict(dict):
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__
    

"""________________________________Before the simulation________________________________"""

@njit(float64(float64))
def VapourPressure(T): # [K]
    T -= 273.15 # [°C]
    return 611.21 * np.exp( (18.678 - T / 234.5) * (T / (257.14 + T)) ) # [Pa]


def InitialCondition(cpar):
    cpar.P_v = VapourPressure(T=cpar.T_inf) # [Pa]
    IC = np.zeros((par.K+3), dtype=np.float64)
    R_0 = cpar.ratio * cpar.R_E
    
    # Equilibrium state
    p_E = cpar.P_inf + 2 * cpar.surfactant * par.sigma / cpar.R_E # [Pa]
    V_E = 4.0 / 3.0 * cpar.R_E**3 * np.pi # [m^3]
    p_gas = p_E - cpar.P_v
    n_gas = p_gas * V_E / (par.R_g * cpar.T_inf) # [mol]
    
    # Isotermic expansion
    V_0 = 4.0 / 3.0 * R_0**3 * np.pi    # [m^3]
    n_H2O = cpar.P_v * V_0 / (par.R_g * cpar.T_inf) # [mol]
    c_H2O = n_H2O / V_0    # [mol/m^3]
    c_gas = n_gas / V_0    # [mol/m^3]

    # Initial conditions
    IC[0] = R_0   # R_0 [m]
    IC[1] = 0.0    # dRdt_0 [m/s]
    IC[2] = cpar.T_inf   # T_0 [K]
    IC[3 + par.indexOfWater] = c_H2O * 1.0e-6    # [mol/cm^3]
    IC[3 + par.indexOfArgon] = c_gas * 1.0e-6    # [mol/cm^3]

    return IC


def Work(cpar):
    cpar.P_v = VapourPressure(T=cpar.T_inf) # [Pa]
    R_0 = cpar.ratio * cpar.R_E # [m]
    V_E = 4.0 / 3.0 * cpar.R_E**3 * np.pi    # [m^3]
    V_0 = 4.0 / 3.0 * R_0**3 * np.pi  # [m^3]
    
    p_E = cpar.P_inf + 2 * cpar.surfactant * par.sigma / cpar.R_E # [Pa]
    p_gas = p_E - cpar.P_v # [Pa]
    n_gas = p_gas * V_E / (par.R_g * cpar.T_inf) # [mol]
    
    W_gas0 = -(cpar.P_v * V_0 + n_gas * par.R_g * cpar.T_inf * np.log(V_0))
    W_gasE = -(cpar.P_v * V_E + n_gas * par.R_g * cpar.T_inf * np.log(V_E))
    W_gas = W_gas0 - W_gasE # [J]
    
    W_surface_tension = par.sigma * cpar.surfactant *  4.0 * np.pi * (R_0**2 - cpar.R_E**2) # [J]
    W_flow = cpar.P_inf * 4.0 / 3.0 * np.pi * (R_0**3 - cpar.R_E**3) # [J]
    return W_gas + W_surface_tension + W_flow    # [J]


"""________________________________Pressures________________________________"""

@njit(float64[:](float64, float64, float64, float64, float64, float64, float64, float64))
def Pressure(R, R_dot, T, T_dot, M, sum_omega_dot, P_inf, surfactant):
    p_inf = P_inf
    p_inf_dot = 0.0
    p = 0.1 * M * par.R_erg * T
    p_dot = p * (sum_omega_dot / M + T_dot / T - 3.0 * R_dot/R)
    p_L = p - (2.0 * surfactant * par.sigma + 4.0 * par.mu_L * R_dot) / R
    p_L_dot = p_dot + (-2.0 * surfactant * par.sigma * R_dot + 4.0 * par.mu_L * R_dot ** 2) / (R ** 2)
    delta = (p_L - p_inf) / par.rho_L
    delta_dot = (p_L_dot - p_inf_dot) / par.rho_L
    return np.array([delta, delta_dot])


"""________________________________NASA polynomials________________________________"""

# Returns the right set of NASA coefficents
@njit(float64[:](int32, float64))
def getCoeffs(k, T):
    a = np.zeros((7), dtype=np.float64)
    if T <= par.TempRange[k][2]: # T <= T_mid
        a = par.a_low[k]
    else:  # T_mid < T
        a = par.a_high[k]
    return a


# returns molar heat capacities, enthalpies and entropies
@njit(float64[:, :](float64[:], float64))
def Thermodynamic(c, T):
    ret = np.zeros((4, par.K), dtype=np.float64)   # [C_p, H, S, C_v]
    exponent = np.array([1, 2, 3, 4, 5])
    for k in range(par.K):
    # get coefficients for T
        a = getCoeffs(k, T)
            
     # calculate
        T_vec = T ** (exponent-1)    # [1, T, T^2, T^3, T^4]
        # Molar heat capacities at constant pressure (isobaric) [erg/mol/K]
        ret[0][k] = par.R_erg * np.sum(a[:par.N] * T_vec)
        # Enthalpies [erg/mol]
        ret[1][k] = par.R_erg * T * (np.sum(a[:par.N] * T_vec / exponent) + a[par.N] / T)
        # Entropies [erg/mol/K]
        ret[2][k] = par.R_erg * (a[0] * np.log(T) + np.sum(a[1:par.N] * T_vec[1:] / (exponent[1:]-1)) + a[par.N+1])
    # Molar heat capacities at constant volume (isochoric) [erg/mol/K]
    ret[3] = ret[0] - par.R_erg 
    return ret


"""________________________________Evaporation________________________________"""

@njit(float64[:](float64, float64, float64, float64, float64, float64))
def Evaporation(M, T, X_H2O, alfa_M, T_inf, P_v):
# condensation and evaporation
    p = 0.1 * M * par.R_erg * T
    p_H2O = X_H2O * p
    n_eva_dot = 1.0e3 * alfa_M * P_v / (par.W[par.indexOfWater] * np.sqrt(2.0 * np.pi * par.R_v * T_inf)) # W_H2O=par.W[5] is in g/mol --> *1000
    n_con_dot = 1.0e3 * alfa_M * p_H2O / (par.W[par.indexOfWater] * np.sqrt(2.0 * np.pi * par.R_v * T))
    n_net_dot = n_eva_dot - n_con_dot
# Molar heat capacity of water at constant volume (isochoric) [erg/mol/K]
    exponent = np.array([1, 2, 3, 4, 5])
    a = getCoeffs(5, T)
    T_vec = T ** (exponent-1)    # [1, T, T^2, T^3, T^4]
    C_v = par.R_erg * (np.sum(a[:par.N] * T_vec) - 1.0)
    
    a = getCoeffs(5, T_inf)
    T_vec = T ** (exponent-1)    # [1, T, T^2, T^3, T^4]
    C_v_inf = par.R_erg * (np.sum(a[:par.N] * T_vec) - 1.0)
# Evaporation energy
    e_eva = C_v_inf * T_inf * 1e-7   # [J/mol]
    e_con = C_v * T * 1e-7
    evap_energy = n_eva_dot * e_eva - n_con_dot * e_con    # [W/m^2]
    return np.array([n_net_dot, evap_energy])


"""________________________________Reaction rates________________________________"""

@njit(float64[:](float64, float64[:]))
def ForwardRate(T, M_eff):
# Reaction rate
    k_forward = par.A * T ** par.b * np.exp(-par.E / (par.R_cal * T))
    
# Pressure dependent reactions
    for j, i in enumerate(par.PressureDependentIndexes):    # i is the number of reaction, j is the index of i's place in par.PressureDependentIndexes
        k_inf = k_forward[i]    # par.A[i] * T ** par.b[i] * np.exp(-par.E[i] / (par.R_cal * T))
        k_0 = par.ReacConst[j][0] * T ** par.ReacConst[j][1] * np.exp(-par.ReacConst[j][2] / (par.R_cal * T))
        P_r = k_0 / k_inf * M_eff[i]
        logP_r = np.log10(P_r)
        
        # Troe mechanis
        F_cent = (1.0 - par.Troe[j][0]) * np.exp(-T / par.Troe[j][1]) + par.Troe[j][0] * np.exp(-T / par.Troe[j][2]) + np.exp(-par.Troe[j][3] / T)
        logF_cent = np.log10(F_cent)
        c = -0.4 - 0.67 * logF_cent
        n = 0.75 - 1.27 * logF_cent
        d = 0.14
        logF = 1.0 / (1.0 + ((logP_r + c) / (n - d * (logP_r + c))) ** 2) * logF_cent
        F = 10.0 ** logF
        
        k_forward[i] = k_inf * P_r / (1.0 + P_r) * F
    return k_forward


@njit(float64[:](float64[:], float64[:], float64[:], float64, float64))
def BackwardRate(k_forward, S, H, T, P_inf):
    DeltaS = np.sum(par.nu * S, axis=1)
    DeltaH = np.sum(par.nu * H, axis=1)
    K_p = np.exp(DeltaS / par.R_erg - DeltaH / (par.R_erg * T))
    K_c = K_p * (P_inf * 10.0 / (par.R_erg * T)) ** np.sum(par.nu, axis=1)
    k_backward = k_forward / K_c
    return k_backward


@njit(float64[:](float64, float64[:], float64[:], float64[:], float64))
def ProductionRate(T, H, S, c, P_inf):
# Third body correction factors
    M_eff = np.sum(c) * np.ones((par.I), dtype = np.float64)    # effective total concentration of the third-body species
    for j, i in enumerate(par.ThirdBodyIndexes):
        M_eff[i] = np.sum(par.alfa[j] * c) 
# Forward and backward rates
    k_forward = ForwardRate(T=T, M_eff=M_eff)
    k_backward = BackwardRate(k_forward=k_forward, S=S, H=H, T=T, P_inf=P_inf)

# Net rates
    q = np.zeros((par.I), dtype = np.float64)
    for i in range(par.I):
        q[i] = k_forward[i] * np.prod(c ** par.nu_forward[i]) - k_backward[i] * np.prod(c ** par.nu_backward[i])
# Third body reactions
    for j, i in enumerate(par.ThirdBodyIndexes):    # i is the number of reaction, j is the index of i in par.ThirdBodyIndexes
        if i not in par.PressureDependentIndexes:
            q[i] *= M_eff[i]
# Production rates
    omega_dot = np.zeros((par.K), dtype=np.float64)
    for k in range(par.K):
        omega_dot[k] = np.sum(par.nu[:, k] * q)
    
    return omega_dot


"""________________________________Differential equation________________________________"""


@njit(float64[:](float64, float64[:], float64, float64, float64, float64, float64))
def f(t, x, alfa_M, T_inf, P_v, P_inf, surfactant):    
    R = x[0]      # bubble radius [m]
    R_dot = x[1]  # [m/s]
    T = x[2]      # temperature [K]
    c = x[3:]     # molar concentration [mol/cm^3]
    M = np.sum(c) # sum of concentration
    X = c / M     # mole fraction [-]
    dxdt = np.zeros(x.shape, dtype = np.float64)
    
# d/dt R
    dxdt[0] = R_dot
# Evaporation
    n_net_dot, evap_energy = Evaporation(M=M, T=T, X_H2O=X[par.indexOfWater], alfa_M=alfa_M, T_inf=T_inf, P_v=P_v)
# Thermodynamics
    [C_p, H, S, C_v] = Thermodynamic(c=c, T=T)
    W_avg = np.sum(X * par.W)
    rho_avg = W_avg * M # or np.sum(c * par.W)
    C_p_avg = np.sum(X * C_p)
    C_v_avg = np.sum(X * C_v)

    lambda_avg = np.sum(X * par.lambdas)
    chi_avg = 10.0 * lambda_avg * W_avg / (C_p_avg * rho_avg)
    l_th = np.inf
    if R_dot != 0.0:
        l_th = np.sqrt(R * chi_avg / abs(R_dot))
    l_th = min(l_th, R / np.pi)
    Q_th_dot = 0.0
    Q_th_dot = lambda_avg * (T_inf - T) / l_th
# d/dt c
    omega_dot = np.zeros((par.K), dtype = np.float64)
    omega_dot = ProductionRate(T=T, H=H, S=S, c=c, P_inf=P_inf)
    c_dot = omega_dot - c * 3.0 * R_dot / R
    c_dot[par.indexOfWater] += 1.0e-6 * n_net_dot * 3.0 / R    # water evaporation
    dxdt[3:] = c_dot
# d/dt T
    Q_r_dot = -np.sum(H * omega_dot)
    p = 0.1 * M * par.R_erg * T
    T_dot = (Q_r_dot + 30.0 / R * (-p * R_dot + Q_th_dot + evap_energy)) / (M * C_v_avg)
    dxdt[2] = T_dot
# d/dt R_dot
    [delta, delta_dot] = Pressure(
        R=R, R_dot=R_dot, T=T, T_dot=T_dot,
        M=M, sum_omega_dot=np.sum(omega_dot),
        P_inf=P_inf, surfactant=surfactant
    )   # delta = (p_L-p_inf) / rho_L
    Nom = (1.0 + R_dot / par.c_L) * delta + R / par.c_L * delta_dot - (1.5 - 0.5 * R_dot / par.c_L) * R_dot ** 2
    Den = (1.0 - R_dot / par.c_L) * R + 4.0 * par.mu_L / (par.c_L * par.rho_L)
    
    dxdt[1] = Nom / Den
    
    return dxdt


"""________________________________Solving________________________________"""

# This funfction solves the differential equation, and returns the numerical solution
# Error codes:
    # 0: succecfully solved with LSODA solver
    # 1: LSODA solver didn't converge, but Radau solver worked
    # 2: LSODA solver had a fatal error, but Radau solver worked
    # 3: Radau solver didn't converge (NO SOLUTION!)
    # 4: Radau solver had a fatal error (NO SOLUTION!)
def solve(cpar, t_int=np.array([0.0, 1.0])):
    error_code = 0
    start = time.process_time()
    cpar.P_v = VapourPressure(T=cpar.T_inf) # [Pa]
    IC = InitialCondition(cpar)
    #f(t, x, alfa_M, T_inf, P_v, P_inf, surfactant)
    
    try:
        num_sol = solve_ivp(f, t_int, IC, method='LSODA', atol = 1e-10, rtol=1e-10, args=(cpar.alfa_M, cpar.T_inf, cpar.P_v, cpar.P_inf, cpar.surfactant))
        if num_sol.success == False: error_code = 1
    except:
        error_code = 2
    if error_code != 0:
        try:
            num_sol = solve_ivp(f, t_int, IC, method='Radau', atol = 1e-10, rtol=1e-10, args=(cpar.alfa_M, cpar.T_inf, cpar.P_v, cpar.P_inf, cpar.surfactant))
            if num_sol.success == False: error = 3
        except:
            error_code = 4
    
    end = time.process_time()
    elapsed_time = (end - start)
    
    if error_code == 4:
        return None, error_code, elapsed_time
    return num_sol, error_code, elapsed_time


"""________________________________Post processing________________________________"""

# This function gets the numerical solution and the control parameters, and returns some datas about the simulation
def getData(cpar, num_sol, error_code, elapsed_time):
    # copy cpar:
    data = dotdict(dict(
        ID=cpar.ID,
        R_E=cpar.R_E,
        ratio=cpar.ratio,
        P_inf=cpar.P_inf,
        alfa_M=cpar.alfa_M,
        T_inf=cpar.T_inf,
        P_v=cpar.P_v,
        surfactant=cpar.surfactant
    ))
    
    # runtime and error
    data.error_code = error_code
    data.elapsed_time = elapsed_time # [s]
    
    # default values
    if error_code==3 or error_code==4:
        data.steps = 0
        loc_min = 0.0
        data.collapse_time = 0.0
        data.T_max = 0.0
        data.x_final = np.zeros((3+par.K), dtype=np.float64)
        data.n_H2 = 0.0
        data.n_O2 = 0.0
        data.work = 0.0
        data.energy = 0.0
    
    # normal functioning
    else:
        data.steps = len(num_sol.t)
        
        # collapse time (first loc min of R)
        loc_min = argrelmin(num_sol.y[:][0])
        data.collapse_time = 0.0
        if not len(loc_min[0]) == 0:
            data.collapse_time = num_sol.t[loc_min[0][0]]
        
        # Energy calculations
        data.T_max = np.max(num_sol.y[:][2]) # maximum of temperature peaks [K]
        data.x_final = num_sol.y[:, -1] # final values of [R, R_dot, T, c_1, ... c_K]
        last_V = 4.0 / 3.0 * (100.0 * data.x_final[0]) ** 3 * np.pi # [cm^3]
        data.n_H2 = data.x_final[3+par.indexOfHydrogen] * last_V # [mol]
        data.n_O2 = data.x_final[3+par.indexOfOxygen] * last_V # [mol]
        m_H2 = 1e-3 * data.n_H2 * par.W[par.indexOfHydrogen] # [kg]
        data.work = Work(cpar) # [J]
        data.energy = 1e-6 * data.work / m_H2 # [MJ/kg]
    
    return data


# keys of data: (except x_final)
keys = ['ID', 'R_E', 'ratio', 'P_inf', 'alfa_M', 'T_inf', 'P_v', 'surfactant', 'error_code', 'elapsed_time', 'steps', 'collapse_time', 'T_max', 'n_H2', 'n_O2', 'work', 'energy']

# This function prints the data dictionary in an organised way
def print_data(data):
    text = f'''Control parameters:
    \tID ={data.ID: .0f}
    \tR_E ={1e6*data.R_E: .2f} [um]
    \tratio ={data.ratio: .2f} [-]
    \tP_inf ={1e-5*data.P_inf: .2f} [bar]
    \talfa_M ={data.alfa_M: .2f} [-]
    \tT_inf ={data.T_inf - 273.15: .2f} [°C]
    \tP_v ={data.P_v: .2f} [Pa]
    \tsurfactant ={data.surfactant: .2f} [bar]
    Simulation info:
    \terror_code ={data.error_code: .0f}
    \telapsed_time ={data.elapsed_time: .2f} [s]
    '''
    
    text += f'''\tsteps ={data.steps: .0f} [-]
    Final state:
    \tR_final ={1e6*data.x_final[0]: .2f} [um];   R_dot_final ={data.x_final[1]} [m/s];   T_final ={data.x_final[2]: .2f} [K]
    \tn_H2 ={data.n_H2} [mol]; n_O2 ={data.n_O2} [mol]
    \tFinal molar concentrations: [mol/cm^3]\n\t\t'''
    
    for k, specie in enumerate(par.species):
        text += f'{specie}: {data.x_final[3+k]};  '
        if (k+1) % 4 == 0: text += f'\n\t\t'
    
    text += f'''\b\bResults:
    \tcollapse_time = {data.collapse_time} [s]
    \tT_max ={data.T_max: .2f} [K]
    \texpansion work = {data.work} [J]
    \thydrogen production ={data.energy: .2f} [MJ/kg]'''
    
    print(text)
    
# This function runs solve() and getData(), then return with data
# input and output is (or can be) normal dictionary
def simulate(cpar, t_int=np.array([0.0, 1.0])):
    cpar = dotdict(cpar)
    num_sol, error_code, elapsed_time = solve(cpar, t_int)
    data = getData(cpar, num_sol, error_code, elapsed_time)
    ret = dict()
    for key in data:
        ret[key] = data[key]
    return ret


"""________________________________Plotting________________________________"""

# This funfction solves the differential equation, and plots it
    # cpar: control parameters in a dictionary
    # t_int: time interval to solve the diffeq in (default: [0, 1] [s])
    #        graphs will be plotted in this intervall, if not default
    # n: how long should the plotted time interval be compared to the collapse time (default: 5 [-])
    # base_name: save plots as .png (default: '' alias do not save)
    #            use base_name='plot' --> plot_1.png, plot_2.png
    #            use base_name='images/plot' to save into images folder
    #            using a folder for images is recommend
    #            this folder have to be created manually
def plot(cpar, t_int=np.array([0.0, 1.0]), n=5.0, base_name=''):
    num_sol, error_code, elapsed_time = solve(cpar, t_int)
    data = getData(cpar, num_sol, error_code, elapsed_time)
    
# Error codes
    if error_code == 0:
        print(f'succecfully solved with LSODA solver')
    if error_code == 1:
        print(f'LSODA solver didn\'t converge, but Radau solver worked')
    if error_code == 2:
        print(f'LSODA solver had a fatal error, but Radau solver worked')
    if error_code == 3:
        print(f'Radau solver didn\'t converge (NO SOLUTION!)')
        return None
    if error_code == 4:
        print(f'Radau solver had a fatal error (NO SOLUTION!)')
        return None
    
# Calculations
    if t_int[1] != 1.0: 
        end_index = -1
    else:
        end_index = np.where(num_sol.t > n * data.collapse_time)[0][0]

    t = num_sol.t[:end_index] # [s]
    R = num_sol.y[0, :end_index] # [m]
    R_dot = num_sol.y[1, :end_index] # [m/s]
    T = num_sol.y[2, :end_index] # [K]
    c = num_sol.y[3:, :end_index] # [mol/cm^3]

    V = 4.0 / 3.0 * (100.0 * R) ** 3 * np.pi # [cm^3]
    n = c * V

# plot R and T
    fig1 = plt.figure(figsize=(20, 6))
    ax1 = fig1.add_subplot(axisbelow=True)
    ax2 = ax1.twinx()
    ax1.plot(t, R / cpar.R_E, color = 'b', linewidth = 1.0)
    ax2.plot(t, T, color = 'r', linewidth = 1.0, linestyle = '-.')

    ax1.set_ylabel('$y1$ [-]')
    ax1.set_xlabel('$t$ [s]')
    ax1.set_ylabel('$R/R_E$ [-]', color = 'b')
    ax2.set_ylabel('$T$ [K]', color = 'r')
    ax1.grid()
    
# textbox with initial conditions
    text = f"""
    Initial conditions:
      {'$R_E$':<25} {1e6*cpar.R_E: .2f}  $[\mu m]$
      {'$R_0/R_E$':<25} {cpar.ratio: .2f}  $[-]$
      {'$P_∞$':<25} {1e-5*cpar.P_inf: .2f}  $[bar]$
      {'$α_M$':<25} {cpar.alfa_M: .2f}  $[-]$
      {'$T_∞$':<25} {cpar.T_inf-273.15: .2f}  $[°C]$
      {'$P_{vapour}$':<25} {cpar.P_v: .1f}  $[Pa]$
      {'$surfactant$':<25} {cpar.surfactant: .2f}  $[-]$
    """
    ax1.text(
        0.75, 0.95, # coordinates
        text, transform=ax1.transAxes,
        horizontalalignment='left', verticalalignment='top',
        fontsize=16, fontstyle='oblique',
        bbox={'facecolor': 'white', 'alpha': 1.0, 'pad': 10},
    )
    
    plt.show()

# plot reactions
    fig2 = plt.figure(figsize=(20, 9))
    ax = fig2.add_subplot(axisbelow=True)

    plt.ylim([1e-23, 1e-11])
    ax.set_yscale('log')
    #O,    H,    H2,     OH,   O2,      H2O,   HO2, H2O2,   O3,  OH_ex
    for i, element in enumerate(par.species):
        ax.plot(t, n[i], label = '$' + element + '$', linewidth = 1.0)

    ax.set_xlabel('$t$ [s]')
    ax.set_ylabel('$n_k$ [mol]')
    ax.grid()
    ax.legend()

    plt.show()
    
# saving the plots
    if base_name != '':
        try:
            fig1.savefig(base_name + '_1.png')
            fig2.savefig(base_name + '_2.png')
        except:
            print(f'Error in saving {base_name}_1.png')

# print data
    print_data(data)
    return None


"""________________________________Save to CSV________________________________"""

class Make_dir:
    # constructor
    def __init__(self, folder_name, file_base_name='output_', seperator=','):
        self.folder_name = folder_name
        self.file_base_name = file_base_name
        self.seperator = seperator
        self.parent_dir = os.getcwd()
        self.save_dir = os.path.join(self.parent_dir, folder_name)
        self.number = 1    # uniqe ID number for csv files
        self.lines = 0    # number of data lines in currently opened CCSV
        self.is_opened = False    # True, if a CSV is opened
        
        if os.path.exists(self.save_dir):
            self.new = False
            self.number = len([1 for file in os.listdir(self.save_dir) if file[-4:] == '.csv']) + 1
            print(f'Folder already exists with {self.number-1} csv in it')
        else:
            self.new = True
            os.mkdir(self.save_dir)
    
    # makes a string from a list e.g. [1, 2, 3] -> '1,2,3'
    def list_to_string(self, array):
        line = ''
        for element in array:
            line += str(element) + self.seperator
        return line[:-1]
        
    # writes a data dict into the currently opened file
    def write_line(self, data):
        line = self.list_to_string([data[key] for key in keys])
        line += self.seperator + self.list_to_string([x for x in data['x_final']])
        self.file.write(line + '\n')
        self.lines += 1
        
    # saves a numerical solution
    def write_solution(self, data, num_sol, file_base_name):
    # create file containing data
        file = os.path.join(self.save_dir, file_base_name + '_data.csv')
        file = open(file, 'w')
        # write header line
        line = self.list_to_string(keys + ['R_last', 'R_dot_last', 'T_last'] + ['c_' + specie + '_last' for specie in par.species])
        file.write(line + '\n')
        # write data
        line = self.list_to_string([data[key] for key in keys])
        line += self.seperator + self.list_to_string([x for x in data.x_final])
        file.write(line + '\n')
        file.close()
    # create file containing num_sol
        file = os.path.join(self.save_dir, file_base_name + '_num_sol.csv')
        file = open(file, 'w')
        # write header line
        line = self.list_to_string(['t', 'R', 'R_dot', 'T'] + ['c_' + specie for specie in par.species])
        file.write(line + '\n')
        # write data
        for i in range(len(num_sol.t)):
            line = self.list_to_string([num_sol.t[i]] + list(num_sol.y[:, i]))
            file.write(line + '\n')
        file.close()
    
    # create new file
    def new_file(self):
        if self.is_opened:
            return None
        file = os.path.join(self.save_dir, self.file_base_name + str(self.number) + '.csv')
        self.file = open(file, 'w')
        self.is_opened = True
        self.number += 1
        self.lines = 0
        # write header line:
        line = self.list_to_string(keys + ['R_last', 'R_dot_last', 'T_last'] + ['c_' + specie + '_last' for specie in par.species])
        self.file.write(line + '\n')
    
    # close file
    def close(self):
        if self.is_opened:
            self.file.close()
            self.is_opened = False
        
    # destructor
    def __del__(self):
        if self.is_opened:
            self.file.close()