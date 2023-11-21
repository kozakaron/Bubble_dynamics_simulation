"""________________________________Libraries________________________________"""

from termcolor import colored
import matplotlib.pyplot as plt   # for plotting
import numpy as np   # matrices, math
from scipy.integrate import solve_ivp   # differential equation solver
from scipy.signal import argrelmin   # loc min finding
import time   # runtime measurement
from datetime import datetime   # for accessing current datetime
import socket   # for accessing computer name
import psutil   # get system information
from numba import njit   # Just In Time compiler
from numba.types import unicode_type, float64, float32, int64, int32   # JIT types
from func_timeout import func_timeout, FunctionTimedOut   # for timeout
import os    # file management
import importlib   # for reloading your own files

# my own files:
try:
    import parameters as par   # numeric constants and coefficents
    importlib.reload(par)   # reload changes you made
except:
    print(print(colored('Error, \'parameters.py\' not found','red')))
try:
    import excitation
    importlib.reload(excitation)
except:
    try:
        import Bubble_dynamics_simulation.excitation as excitation
        importlib.reload(excitation)
    except:
        print(colored(f'Error, \'excitation.py\' not found', 'red'))    

"""________________________________Settings________________________________"""

enable_heat_transfer = True
enable_evaporation = False
enable_reactions = True
enable_dissipated_energy = False
target_specie = 'NH3' # Specie to calculate energy effiqiency
excitation_type = 'no_excitation' # function to calculate pressure excitation

"""________________________________General________________________________"""

class dotdict(dict):
    """Dot notation access to dictionary attributes. 
    Instead of dictionary['key'] you can use dictionary.key"""
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

def copy(input):
    """deep copy of input"""

    if type(input) == list:
        return [copy(element) for element in input]
    elif type(input) == dict:
        return {key: copy(value) for key, value in input.items()}
    elif type(input) == np.ndarray:
        return input.copy()
    elif type(input) == dotdict:
        return dotdict({key: copy(value) for key, value in input.items()})
    else:
        return input

"""________________________________Before the simulation________________________________"""

def colorTF(boolean):
    return colored(str(boolean), 'green') if boolean else colored(str(boolean), 'red')

Excitation, excitation_args, excitation_units = excitation.getExcitation(excitation_type=excitation_type)
if par.indexOfWater == -1:
    enable_evaporation = False
print(f'model: {par.model}')
print(f'target specie: {target_specie}')
print(f'excitation: {excitation_type} (control parameters: {excitation_args})')
print(f'enable heat transfer: {colorTF(enable_heat_transfer)}\tenable evaporation: {colorTF(enable_evaporation)}\tenable reactions: {colorTF(enable_reactions)}\tenable dissipated energy: {colorTF(enable_dissipated_energy)}')
if target_specie not in par.species:
    print(colored(f'Error, target specie \'{target_specie}\' not found in parameters.py', 'red'))
    
def example_cpar(normal_dict=False):
    '''Provides an example of the control parameter dictionary. Use print_cpar() to print it. Parameters:
    * normal_dict: if True, returns a normal dictionary, else returns a dotdict
    
    Returns:
    * cpar: control parameter dictionary'''
    
    cpar = dict(
        ID = 0,                            # ID of control parameter (not used during calculation)
    # Initial conditions:
        R_E = 10.0e-6,                     # bubble equilibrium radius [m]
        ratio = 1.0,                       # initial radius / equilibrium radius R_0/R_E [-]
        gases = [0],                       # indexes of species in initial bubble (list of species indexes)
        fractions = [1.0],                 # molar fractions of species in initial bubble (list of fractions for every gas)
    # Ambient parameters:
        P_amb = 1.0 * par.atm2Pa,          # ambient pressure [Pa]
        T_inf = 20.0 + par.absolute_zero,  # ambient temperature [K]
    # Liquid parameters:
        alfa_M = par.alfa_M,               # water accommodation coefficient [-]
        P_v = par.P_v,                     # vapour pressure [Pa]
        mu_L = par.mu_L,                   # dynamic viscosity [Pa*s]
        c_L = par.c_L,                     # sound speed [m/s]
        surfactant = 1.00,                 # surfactant (surface tension modfier) [-]
    )

    for arg in excitation_args:
        cpar[arg] = 0.0
    if target_specie == 'NH3':
        cpar['gases'] = [par.index['H2'], par.index['N2']]
        cpar['fractions'] = [0.75, 0.25]
    else:
        cpar['gases'] = [par.index['O2']]
        cpar['fractions'] = [1.0]

    if normal_dict:
        return cpar
    else:
        return dotdict(cpar)

@njit(float64(float64))
def VapourPressure(T): # [K]
    T -= 273.15 # [°C]
    return 611.21 * np.exp( (18.678 - T / 234.5) * (T / (257.14 + T)) ) # [Pa]

@njit(float64(float64))
def Viscosity(T): # [K], pressure dependence is neglected
    return 1.856e-14 * np.exp(4209.0/T + 0.04527*T - 3.376e-5*T**2) # [Pa*s]

def InitialCondition(cpar, evaporation=False):
    if not 'P_v' in cpar:
        cpar.P_v = VapourPressure(T=cpar.T_inf) # [Pa]
    if not 'mu_L' in cpar:
        cpar.mu_L = Viscosity(T=cpar.T_inf) # [Pa]
    if not 'c_L' in cpar:
        cpar.c_L = par.c_L # [m/s]
    if type(cpar.gases) != list:
        cpar.gases = [cpar.gases]
    if type(cpar.fractions) != list:
        cpar.fractions = [cpar.fractions]
    if round(sum(cpar.fractions), 5) != 1.0:
        print(print(colored(f'Warning, in InitialCondition(), sum of cpar.fractions isn\'t 1: {cpar.fractions}','yellow')))
    if len(cpar.gases) != len(cpar.fractions):
        print(print(colored(f'Warning, in InitialCondition(), len(cpar.gases) != len(cpar.fractions): {cpar.gases} != {cpar.fractions}','yellow')))
    IC = np.zeros((par.K+4), dtype=np.float64)
    R_0 = cpar.ratio * cpar.R_E
    
    # Equilibrium state
    p_E = cpar.P_amb + 2.0 * cpar.surfactant * par.sigma / cpar.R_E # [Pa]
    V_E = 4.0 / 3.0 * cpar.R_E**3 * np.pi # [m^3]
    p_gas = p_E - cpar.P_v if evaporation else p_E
    lowpressure_error = lowpressure_warning = False
    if p_gas < 0.0:
        #print(colored('Error! The pressure of the gas is negative!', 'red'))
        lowpressure_error = True
    n_gas = p_gas * V_E / (par.R_g * cpar.T_inf) # [mol]
    
    # Isotermic expansion
    V_0 = 4.0 / 3.0 * R_0**3 * np.pi    # [m^3]
    n_H2O = cpar.P_v * V_0 / (par.R_g * cpar.T_inf) if evaporation else 0.0 # [mol]
    c_H2O = n_H2O / V_0    # [mol/m^3]
    c_gas = n_gas / V_0    # [mol/m^3]
    p_gas = c_gas * par.R_g * cpar.T_inf # [Pa]
    P_amb_min = cpar.P_v if evaporation else 0.0 # [Pa]
    P_amb_min += p_gas - 2.0 * cpar.surfactant * par.sigma / R_0 # [Pa]
    if P_amb_min < cpar.P_v:
        #print(colored('Warning! The pressure during the expansion is lower, than the saturated water pressure', 'yellow'))
        lowpressure_warning = True

    # Initial conditions
    IC[0] = R_0   # R_0 [m]
    IC[1] = 0.0    # dRdt_0 [m/s]
    IC[2] = cpar.T_inf   # T_0 [K]
    if evaporation and cpar.indexOfWater != -1:
        IC[3 + par.index['H2O']] = c_H2O * 1.0e-6    # [mol/cm^3]
    for index, fraction in zip(cpar.gases, cpar.fractions):
        IC[3 + index] = fraction * c_gas * 1.0e-6    # [mol/cm^3]
    IC[3 + par.K] = 0.0 # dissipated acoustic energy [J]

    return IC, lowpressure_error, lowpressure_warning


def Work(cpar, evaporation=False):
    if not 'P_v' in cpar:
        cpar.P_v = VapourPressure(T=cpar.T_inf) # [Pa]
    if not 'mu_L' in cpar:
        cpar.mu_L = Viscosity(T=cpar.T_inf) # [Pa]
    R_0 = cpar.ratio * cpar.R_E # [m]
    V_E = 4.0 / 3.0 * cpar.R_E**3 * np.pi    # [m^3]
    V_0 = 4.0 / 3.0 * R_0**3 * np.pi  # [m^3]
    
    p_E = cpar.P_amb + 2 * cpar.surfactant * par.sigma / cpar.R_E # [Pa]
    p_gas = p_E - cpar.P_v if evaporation else p_E # [Pa]
    n_gas = p_gas * V_E / (par.R_g * cpar.T_inf) # [mol]
    
    W_gas0 = -(cpar.P_v * V_0 + n_gas * par.R_g * cpar.T_inf * np.log(V_0)) if evaporation else -(n_gas * par.R_g * cpar.T_inf * np.log(V_0))
    W_gasE = -(cpar.P_v * V_E + n_gas * par.R_g * cpar.T_inf * np.log(V_E)) if evaporation else -(n_gas * par.R_g * cpar.T_inf * np.log(V_E))
    W_gas = W_gas0 - W_gasE # [J]
    
    W_surface_tension = par.sigma * cpar.surfactant *  4.0 * np.pi * (R_0**2 - cpar.R_E**2) # [J]
    W_flow = cpar.P_amb * 4.0 / 3.0 * np.pi * (R_0**3.0 - cpar.R_E**3.0) # [J]
    return W_gas + W_surface_tension + W_flow    # [J]


"""________________________________Pressures________________________________"""

@njit(float64[:](float64, float64, float64, float64, float64, float64, float64, float64, float64[:]))
def Pressure(t, R, R_dot, mu_L, surfactant, p, p_dot, P_amb, args):
    p_Inf, p_Inf_dot = Excitation(t, P_amb, args)
    p_L = p - (2.0 * surfactant * par.sigma + 4.0 * mu_L * R_dot) / R
    p_L_dot = p_dot + (-2.0 * surfactant * par.sigma * R_dot + 4.0 * mu_L * R_dot ** 2) / (R ** 2)
    delta = (p_L - p_Inf) / par.rho_L
    delta_dot = (p_L_dot - p_Inf_dot) / par.rho_L
    return np.array([delta, delta_dot], dtype=np.float64)


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
        ret[1][k] = par.R_erg * ( T * np.sum(a[:par.N] * T_vec / exponent) + a[par.N])
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
    n_eva_dot = 1.0e3 * alfa_M * P_v / (par.W[par.indexOfWater] * np.sqrt(2.0 * np.pi * par.R_v * T_inf))
    n_con_dot = 1.0e3 * alfa_M * p_H2O / (par.W[par.indexOfWater] * np.sqrt(2.0 * np.pi * par.R_v * T))
    n_net_dot = n_eva_dot - n_con_dot
# Molar heat capacity of water at constant volume (isochoric) [erg/mol/K]
    exponent = np.array([1, 2, 3, 4, 5])
    a = getCoeffs(par.indexOfWater, T)
    T_vec = T ** (exponent-1)    # [1, T, T^2, T^3, T^4]
    C_v = par.R_erg * (np.sum(a[:par.N] * T_vec) - 1.0)
    
    a = getCoeffs(par.indexOfWater, T_inf)
    T_vec = T_inf ** (exponent-1)    # [1, T, T^2, T^3, T^4]
    C_v_inf = par.R_erg * (np.sum(a[:par.N] * T_vec) - 1.0)
# Evaporation energy
    e_eva = C_v_inf * T_inf * 1e-7   # [J/mol]
    e_con = C_v * T * 1e-7
    evap_energy = n_eva_dot * e_eva - n_con_dot * e_con    # [W/m^2]
    
    return np.array([n_net_dot, evap_energy], dtype=np.float64)


"""________________________________Reaction rates________________________________"""

@njit(float64[:](float64, float64[:], float64))
def ForwardRate(T, M_eff, M):
# Reaction rate
    k_forward = par.A * T ** par.b * np.exp(-par.E / (par.R_cal * T))
    
# Pressure dependent reactions
    for j, i in enumerate(par.PressureDependentIndexes):    # i is the number of reaction, j is the index of i's place in par.PressureDependentIndexes
        k_inf = k_forward[i]    # par.A[i] * T ** par.b[i] * np.exp(-par.E[i] / (par.R_cal * T))
        k_0 = par.ReacConst[j][0] * T ** par.ReacConst[j][1] * np.exp(-par.ReacConst[j][2] / (par.R_cal * T))
        P_r = k_0 / k_inf * M_eff[i]
        logP_r = np.log10(P_r)
        
         # Lindemann formalism
        if i in par.LindemannIndexes:
            F = 1.0

        # Troe formalism
        elif i in par.TroeIndexes:
            F_cent = (1.0 - par.Troe[j][0]) * np.exp(-T / par.Troe[j][1]) + par.Troe[j][0] * np.exp(-T / par.Troe[j][2]) + np.exp(-par.Troe[j][3] / T)
            logF_cent = np.log10(F_cent)
            c2 = -0.4 - 0.67 * logF_cent
            n = 0.75 - 1.27 * logF_cent
            d = 0.14
            logP_r = np.log10(P_r)
            logF = 1.0 / (1.0 + ((logP_r + c2) / (n - d * (logP_r + c2))) ** 2) * logF_cent
            F = 10.0 ** logF
        
        # SRI formalism
        elif i in par.SRIIndexes: 
            X = 1.0 / (1.0 + np.log10(P_r)**2)
            F = par.SRI[j][3] * (par.SRI[j][0] * np.exp(-par.SRI[j][1] / T) + np.exp(-T / par.SRI[j][2]))**X * T ** par.SRI[j][4]
        else:
            print('Error, the pressure-dependent reaction cannot be groupped in any type of pressure-dependent reactions!')
      # Pressure dependent reactions END
    
        k_forward[i] = k_inf * P_r / (1.0 + P_r) * F

  # PLOG reactions
    if par.PlogCount > 0:
        p = 0.1 * M * par.R_erg * T
    for j, i in enumerate(par.PlogIndexes):
        if p < par.Plog[3*j+1][0]:
            k_1 = par.Plog[3*j][1] * T ** par.Plog[3*j][2] * np.exp(-par.Plog[3*j][3] / (par.R_cal * T))
            k_2 = par.Plog[3*j+1][1] * T ** par.Plog[3*j+1][2] * np.exp(-par.Plog[3*j+1][3] / (par.R_cal * T))
            ln_k = np.log(k_1) + (np.log(p) - np.log(par.Plog[3*j][0])) / (np.log(par.Plog[3*j+1][0]) - np.log(par.Plog[3*j][0])) * (np.log(k_2) - np.log(k_1))
            k_forward[i] = np.exp(ln_k)
        else:
            k_2 = par.Plog[3*j+1][1] * T ** par.Plog[3*j+1][2] * np.exp(-par.Plog[3*j+1][3] / (par.R_cal * T))
            k_3 = par.Plog[3*j+2][1] * T ** par.Plog[3*j+2][2] * np.exp(-par.Plog[3*j+2][3] / (par.R_cal * T))
            ln_k = np.log(k_2) + (np.log(p) - np.log(par.Plog[3*j+1][0])) / (np.log(par.Plog[3*j+2][0]) - np.log(par.Plog[3*j+1][0])) * (np.log(k_3) - np.log(k_2))
            k_forward[i] = np.exp(ln_k)
    return k_forward


@njit(float64[:](float64[:], float64[:], float64[:], float64, float64))
def BackwardRate(k_forward, S, H, T, P_amb):
    DeltaS = np.sum(par.nu * S, axis=1)
    DeltaH = np.sum(par.nu * H, axis=1)
    K_p = np.exp(DeltaS / par.R_erg - DeltaH / (par.R_erg * T))
    K_c = K_p * (P_amb * 10.0 / (par.R_erg * T)) ** np.sum(par.nu, axis=1)
    k_backward = k_forward / K_c
    for i in par.IrreversibleIndexes:
        k_backward[i] = 0.0
    return k_backward


@njit(float64[:](float64, float64[:], float64[:], float64[:], float64, float64))
def ProductionRate(T, H, S, c, P_amb, M):
# Third body correction factors
    M_eff = np.sum(c) * np.ones((par.I), dtype = np.float64)    # effective total concentration of the third-body 
    for j, i in enumerate(par.ThirdBodyIndexes):
        M_eff[i] = np.sum(par.alfa[j] * c) 
# Forward and backward rates
    k_forward = ForwardRate(T=T, M_eff=M_eff, M=M)
    k_backward = BackwardRate(k_forward=k_forward, S=S, H=H, T=T, P_amb=P_amb)

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

@njit(float64[:](float64, float64[:], float64, float64, float64, float64, float64, float64, float64, float64[:]))
def f(t, x, P_amb, alfa_M, T_inf, surfactant, P_v, mu_L, c_L, ex_args):   
    R = x[0]      # bubble radius [m]
    R_dot = x[1]  # [m/s]
    T = x[2]      # temperature [K]
    c = x[3:-1]     # molar concentration [mol/cm^3]
    M = np.sum(c) # sum of concentration
    X = c / M     # mole fraction [-]
    dxdt = np.zeros(x.shape, dtype = np.float64)
    
# d/dt R
    dxdt[0] = R_dot
# Evaporation
    n_net_dot = 0.0  
    evap_energy = 0.0
    if enable_evaporation:
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
    if enable_heat_transfer:
        Q_th_dot = lambda_avg * (T_inf - T) / l_th
# d/dt c
    omega_dot = np.zeros((par.K), dtype = np.float64)
    if enable_reactions:
        omega_dot = ProductionRate(T=T, H=H, S=S, c=c, P_amb=P_amb, M=M)
    c_dot = omega_dot - c * 3.0 * R_dot / R
    c_dot[par.indexOfWater] += 1.0e-6 * n_net_dot * 3.0 / R    # water evaporation
    dxdt[3:-1] = c_dot
# d/dt T
    sum_omega_dot = np.sum(omega_dot)
    Q_r_dot = -np.sum(H * omega_dot) + sum_omega_dot * par.R_erg * T
    p = 0.1 * M * par.R_erg * T # Partial pressure of the gases [Pa]
    T_dot = (Q_r_dot + 30.0 / R * (-p * R_dot + Q_th_dot + evap_energy)) / (M * C_v_avg)
    p_dot = p * (sum_omega_dot / M + T_dot / T - 3.0 * R_dot / R) # for later use
    dxdt[2] = T_dot
# d/dt R_dot
    [delta, delta_dot] = Pressure(t=t,
        R=R, R_dot=R_dot, mu_L=mu_L, surfactant=surfactant,
        p=p, p_dot=p_dot, P_amb=P_amb, args=ex_args
    )   # delta = (p_L-P_amb) / rho_L
    
    Nom = (1.0 + R_dot / c_L) * delta + R / c_L * delta_dot - (1.5 - 0.5 * R_dot / c_L) * R_dot ** 2
    Den = (1.0 - R_dot / c_L) * R + 4.0 * mu_L / (c_L * par.rho_L)
    
    dxdt[1] = Nom / Den
    
    if enable_dissipated_energy:
        V_dot=4.0 * R * R * R_dot * np.pi
        integrand_th = -(p * (1 + R_dot / c_L) + R / c_L * p_dot) * V_dot
        integrand_v = 16.0 * np.pi * mu_L * (R * R_dot*R_dot + R * R * R_dot * dxdt[1] / c_L)
        integrand_r = 4.0 * np.pi / c_L * R * R * R_dot * (R_dot * p + p_dot * R - 0.5 * par.rho_L * R_dot * R_dot * R_dot - par.rho_L * R * R_dot * dxdt[1])

        dxdt[-1]=(integrand_th + integrand_v + integrand_r)
    else:
        dxdt[-1]=0.0
    
    return dxdt


"""________________________________Solving________________________________"""

def solve(cpar, t_int=np.array([0.0, 1.0]), LSODA_timeout=30, Radau_timeout=300):
    """
    This funfction solves the differential equation, and returns the numerical solution.
    Parameters:
     * cpar: control parameters
     * t_int: time interval
     * LSODA_timeout: timeout for LSODA solver
     * Radau_timeout: timeout for Radau solver

    Returns:
     * num_sol: numerical solution. use num_sol.t and num_sol.y to get the time and the solution
     * error_code: see de.error_codes: dict, de.get_errors()
     * elapsed_time: elapsed time
    """
    error_code = 0
    start = time.time()
    if not 'P_v' in cpar:
        cpar.P_v = VapourPressure(T=cpar.T_inf) # [Pa]
    if not 'mu_L' in cpar:
        cpar.mu_L = Viscosity(T=cpar.T_inf) # [Pa]
    ex_args = []
    for name, unit in zip(excitation_args, excitation_units):
        if name in cpar:
            ex_args.append(cpar[name])
        else:
            ex_args.append(0.0)
            print(colored(f'Warning! Pressure excitation argument \'{name} [{unit}]\' is not in cpar. 0 is used instead. ', 'yellow'))
    ex_args = np.array(ex_args, dtype=np.float64)
    IC, lowpressure_error, lowpressure_warning = InitialCondition(cpar, enable_evaporation)
    if lowpressure_error:
        error_code += 100
        return None, error_code, 0.0
    elif lowpressure_warning:
        error_code += 200
    
    # solving d/dt x=f(t, x, cpar)
    try: # try-catch block
        num_sol = func_timeout( # timeout block
            LSODA_timeout, solve_ivp,
            kwargs=dict(fun=f, t_span=t_int, y0=IC, method='LSODA', atol = 1e-10, rtol=1e-10, # solve_ivp's arguments
                        args=(cpar.P_amb, cpar.alfa_M, cpar.T_inf, cpar.surfactant, cpar.P_v, cpar.mu_L, cpar.c_L, ex_args) # f's arguments
            )
        )
        if num_sol.success == False:
            error_code += 1
    except FunctionTimedOut:
        error_code += 2
    except:
        error_code += 3
    if error_code % 10 != 0:
        try: # try-catch block
            num_sol = func_timeout( # timeout block
                Radau_timeout, solve_ivp, 
                kwargs=dict(fun=f, t_span=t_int, y0=IC, method='Radau', atol = 1e-10, rtol=1e-10, # solve_ivp's arguments
                            args=(cpar.P_amb, cpar.alfa_M, cpar.T_inf, cpar.surfactant, cpar.P_v, cpar.mu_L, cpar.c_L, ex_args)) # f's arguments
            )
            if num_sol.success == False:
                error_code += 40
        except FunctionTimedOut:
            error_code += 50
        except:
            error_code += 60
    
    end = time.time()
    elapsed_time = (end - start)
    
    if error_code % 100 > 3:
        return None, error_code, elapsed_time
    return num_sol, error_code, elapsed_time

# error codes description
error_codes = { # this is also a dictionary
    'xx0': dict(describtion='succecfully solved with LSODA solver', color='green'),
    'xx1': dict(describtion='LSODA solver didn\'t converge', color='yellow'),
    'xx2': dict(describtion='LSODA solver timed out', color='yellow'),
    'xx3': dict(describtion='LSODA solver had a fatal error', color='yellow'),
    'x0x': dict(describtion='succecfully solved with Radau solver', color='green'),
    'x4x': dict(describtion='Radau solver didn\'t converge (NO SOLUTION!)', color='red'),
    'x5x': dict(describtion='Radau solver timed out (NO SOLUTION!)', color='red'),
    'x6x': dict(describtion='Radau solver had a fatal error (NO SOLUTION!)', color='red'),
    '1xx': dict(describtion='Low pressure error: The pressure of the gas is negative', color='red'),
    '2xx': dict(describtion='Low pressure warning: The pressure during the expansion is lower, than the saturated water pressure', color='yellow'),
}

def get_errors(error_code, printit=False):
    '''
    * Input: error_code (int) 
    * Output: list of error codes (str)
    * Also prints colored errors, if printit=True
    '''
    # get digits of error_code
    first_digit = f'xx{error_code % 10}'
    second_digit = f'x{(error_code // 10) % 10}x'
    third_digit = f'{(error_code // 100) % 10}xx'

    # get only relevant errors
    errors = []
    errors.append(first_digit)
    if first_digit != 'xx0':
        errors.append(second_digit)
    if third_digit != '0xx':
        errors.append(third_digit)

    # determine if soultion was succesfull
    success = '1xx' not in errors and ('xx0' in errors or 'x0x' in errors)

    # print errors
    if printit:
        for error in errors:
            print(colored(error_codes[error]['describtion'], error_codes[error]['color']))
    return errors, success

"""________________________________Post processing________________________________"""

# This function gets the numerical solution and the control parameters, and returns some datas about the simulation
def get_data(cpar, num_sol, error_code, elapsed_time):
    # copy cpar:
    data = dotdict(dict(
        ID=cpar.ID,
        R_E=cpar.R_E,
        ratio=cpar.ratio,
        P_amb=cpar.P_amb,
        alfa_M=cpar.alfa_M,
        T_inf=cpar.T_inf,
        P_v=cpar.P_v,
        mu_L=cpar.mu_L,
        surfactant=cpar.surfactant,
        gases=cpar.gases,
        fractions=cpar.fractions,
        c_L=cpar.c_L,
    ))
    for name, unit in zip(excitation_args, excitation_units):
        if name in cpar:
            data[name] = cpar[name]
        else:
            data[name] = 0.0
            #print(colored(f'Warning! Pressure excitation argument \'{name} [{unit}]\' is not in cpar. 0 is used instead. ', 'yellow'))
    
    # runtime and error
    data.error_code = error_code
    data.elapsed_time = elapsed_time # [s]
    
    # default values
    data.steps = 0
    data.collapse_time = 0.0
    data.T_max = 0.0
    data.x_initial = np.zeros((4+par.K), dtype=np.float64)
    data.x_final = np.zeros((4+par.K), dtype=np.float64)
    data[f'n_{target_specie}'] = 0.0
    data.m_target = 0.0
    data.expansion_work = 0.0
    data.dissipated_acoustic_energy = 0.0
    data.energy_efficiency = 1.0e30
    data.enable_heat_transfer = enable_heat_transfer
    data.enable_evaporation = enable_evaporation
    data.enable_reactions = enable_reactions
    data.enable_dissipated_energy = enable_dissipated_energy
    data.excitation_type = excitation_type
    data.target_specie = target_specie
    errors, success = get_errors(error_code)
    data.success = success
    if not success:
        return data
    
    # normal functioning
    data.steps = len(num_sol.t)
    data.x_initial = num_sol.y[:, 0] # initial values of [R, R_dot, T, c_1, ... c_K]
        
    # collapse time (first loc min of R)    TODO fix
    loc_min = argrelmin(num_sol.y[:][0])
    data.collapse_time = 0.0
    if not len(loc_min[0]) == 0:
        data.collapse_time = num_sol.t[loc_min[0][0]]
        
    # Energy calculations
    data.T_max = np.max(num_sol.y[:][2]) # maximum of temperature peaks [K]
    data.x_final = num_sol.y[:, -1] # final values of [R, R_dot, T, c_1, ... c_K]
    last_V = 4.0 / 3.0 * (100.0 * data.x_final[0]) ** 3 * np.pi # [cm^3]
    data[f'n_{target_specie}'] = data.x_final[3+par.index[target_specie]] * last_V # [mol]
    m_target = 1.0e-3 * data[f'n_{target_specie}'] * par.W[par.index[target_specie]] # [kg]
    data.expansion_work = Work(cpar, enable_evaporation) # [J]
    data.dissipated_acoustic_energy = data.x_final[-1]  # [J]
    all_work = data.expansion_work + data.dissipated_acoustic_energy
    data.energy_efficiency = 1.0e-6 * all_work / m_target if m_target > 0.0 else 1.0e30 # [MJ/kg]
    data.target_specie = target_specie
    return data

# keys of data: (except x_final)
keys = ['ID', 'R_E', 'ratio', 'P_amb', 'alfa_M', 'T_inf', 'P_v', 'mu_L', 'gases', 'fractions', 'surfactant', 'c_L',
        'error_code', 'success', 'elapsed_time', 'steps', 'collapse_time', 'T_max', f'n_{target_specie}', 'expansion_work', 'dissipated_acoustic_energy', 'energy_efficiency',
        'enable_heat_transfer', 'enable_evaporation', 'enable_reactions', 'enable_dissipated_energy', 'excitation_type', 'target_specie'] + excitation_args

# used in print_cpar
def print_line(name, value, comment, print_it=False):
    text = f'    {name} = '
    if type(value) == str:
        text += f'\'{value}\','
    elif type(value) == int:
        text += f'{value},'
    elif type(value) == float:
        if abs(value) < 1e-12:
            text += f'0.0,'
        elif abs(value) < 1e-6:
            text += f'{value: .6e},'
        elif abs(value) < 1e-3:
            text += f'{value: .8f},'
        elif abs(value) < 1.0:
            text += f'{value: .4f},'
        else:
            text += f'{value: .2f},'
    elif type(value) == np.ndarray or type(value) == list:
        if name == 'gases':
            gases= ''.join([f'par.index[\'{par.species[i]}\'], ' for i in value])
            text += f'[{gases[:-2]}],'
        if name == 'fractions':
            fractions = ''.join([f'{i}, ' for i in value])
            text += f'[{fractions[:-2]}],'
    else:
        text += f'{value},'

    if print_it:
        print(f'{text: <48} # {comment}')
    else:
        return f'{text: <48} # {comment}\n'

# print cpar in an organised way (works with dict, dotdict)
def print_cpar(cpar, without_code=False, print_it=True):
    text = ''
    if not without_code:
        text += f'cpar = de.dotdict(dict(\n'
    text += print_line('ID', int(cpar['ID']), 'ID of control parameter (not used during calculation)')
    text += f'  # Initial conditions:\n'
    text += print_line('R_E', float(cpar['R_E']), 'bubble equilibrium radius [m]')
    text += print_line('ratio', float(cpar['ratio']), 'initial radius / equilibrium radius R_0/R_E [-]')
    text += print_line('gases', cpar['gases'], 'indexes of species in initial bubble (list of species indexes)')
    text += print_line('fractions', cpar['fractions'], 'molar fractions of species in initial bubble (list of fractions for every gas)')
    text += f'  # Ambient parameters:\n'
    text += print_line('P_amb', float(cpar['P_amb']), 'ambient pressure [Pa]')
    text += print_line('T_inf', float(cpar['T_inf']), 'ambient temperature [K]')
    text += f'  # Liquid parameters:\n'
    text += print_line('alfa_M', float(cpar['alfa_M']), 'water accommodation coefficient [-]')
    text += print_line('P_v', float(cpar['P_v']), 'vapour pressure [Pa]')
    text += print_line('mu_L', float(cpar['mu_L']), 'dynamic viscosity [Pa*s]')
    text += print_line('c_L', float(cpar['c_L']), 'sound speed [m/s]')
    text += print_line('surfactant', float(cpar['surfactant']), 'surfactant (surface tension modfier) [-]')
    text += f'  # Excitation parameters: (excitation_type = {excitation_type})\n'
    for arg, unit in zip(excitation_args, excitation_units):
        text += print_line(arg, cpar[arg], f'[{unit}]')
    if not without_code:
        text += f'))\n\n# Calculate pressure/temperature dependent parameters:\n'
        text += f'cpar.mu_L = de.Viscosity(cpar.T_inf)\n'
        text += f'cpar.P_v = de.VapourPressure(cpar.T_inf)\n'
    if print_it:
        print(text)
    else:
        return text

# This function prints the data dictionary in an organised way
def print_data(cpar, data, print_it=True):
    text = f'Control parameters:\n'
    text += print_cpar(cpar, without_code=True, print_it=False)
    text += f'''\nSimulation info:
    error_code ={data.error_code: .0f} (success = {data.success})
    elapsed_time ={data.elapsed_time: .2f} [s]
    steps ={data.steps: .0f} [-]'''
    
    text += f'''\nFinal state:
    R_final ={1e6*data.x_final[0]: .2f} [um];   R_dot_final ={data.x_final[1]} [m/s];   T_final ={data.x_final[2]: .2f} [K]
    n_{target_specie}_final ={data[f'n_{target_specie}']: .2e} [mol]
    Final molar concentrations: [mol/cm^3]\n        '''
    
    for k, specie in enumerate(par.species):
        text += f'{specie: <6}: {data.x_final[3+k]: 24};    '
        if (k+1) % 4 == 0: text += f'\n        '
    
    text += f'''\nResults:
    collapse_time = {data.collapse_time} [s]
    T_max ={data.T_max: .2f} [K]
    expansion work = {data.expansion_work} [J]
    dissipated acoustic energy = {data.dissipated_acoustic_energy} [J]
    energy efficiency = {data.energy_efficiency} [MJ/kg of {target_specie}]'''
    
    if print_it:
        print(text)
    else:
        return text
    
# This function runs solve() and get_data(), then return with data
# input and output is (or can be) normal dictionary
# it is used for multithreading (e.g. in bruteforce_parameter_study.inp)
def simulate(kwargs):
    args = dict(t_int=np.array([0.0, 1.0]), LSODA_timeout=30, Radau_timeout=300)
    for key in kwargs:
        args[key] = kwargs[key]
    args = dotdict(args)
    cpar = dotdict(args.cpar)
    num_sol, error_code, elapsed_time = solve(cpar, args.t_int, LSODA_timeout=args.LSODA_timeout, Radau_timeout=args.Radau_timeout)
    data = get_data(cpar, num_sol, error_code, elapsed_time)
    return dict(data)


"""________________________________Plotting________________________________"""

def plot(cpar, t_int=np.array([0.0, 1.0]), n=5.0, base_name='', LSODA_timeout=30, Radau_timeout=300, presentation_mode=False, plot_pressure=False, show_legend=False, show_cpar=True):
    """
    This funfction solves the differential equation, and plots it.
    Parameters:
     * cpar: control parameters in a dictionary
     * t_int: time interval to solve the diffeq in (default: [0, 1] [s])
           graphs will be plotted in this intervall, if not default
     * n: how long should the plotted time interval be compared to the collapse time (default: 5 [-])
     * base_name: save plots as .png (default: '' alias do not save)
               use base_name='plot' --> plot_1.png, plot_2.png
               use base_name='images/plot' to save into images folder
               using a folder for images is recommend
               this folder have to be created manually
     * LSODA_timeout, Radau_timeout: timeout (maximum runtime) for different solvers in solve() in seconds
     * presentation_mode: if True, the plot will be in presentation mode (default: False)
     * plot_pressure: if True, the pressure will be plotted (default: False)
     * show_legend: if True, the legend will be visible with every single species (default: False)
     * show_cpar: if True, the control parameters will be printed on the plot (default: False)
    """
    if type(cpar) == dict:
        cpar = dotdict(cpar)

    num_sol, error_code, elapsed_time = solve(cpar, t_int, LSODA_timeout, Radau_timeout)
    data = get_data(cpar, num_sol, error_code, elapsed_time)
    
# Print errors
    errors, success = get_errors(error_code, printit=True)
    if not success:
        print_data(cpar, data)
        return None
    
# Calculations
    if t_int[1] != 1.0: 
        end_index = -1
    else:
        end_index = np.where(num_sol.t > n * data.collapse_time)[0][0]

    if num_sol.t[end_index] < 1e-3:
        t = num_sol.t[:end_index] * 1e6 # [us]
    else:
        t = num_sol.t[:end_index] * 1e3 # [ms]
    R = num_sol.y[0, :end_index] # [m]
    R_dot = num_sol.y[1, :end_index] # [m/s]
    T = num_sol.y[2, :end_index] # [K]
    c = num_sol.y[3:, :end_index] # [mol/cm^3]

    V = 4.0 / 3.0 * (100.0 * R) ** 3 * np.pi # [cm^3]
    n = c * V
    if plot_pressure:
        internal_pressure = 1e-6 * cpar.P_v + np.sum(n, axis=0) * par.R_g * T / V # [MPa]

# plot R and T
    linewidth = 2.0 if presentation_mode else 1.0
    plt.rcParams.update({'font.size': 24 if presentation_mode else 18})
    fig1 = plt.figure(figsize=(16, 9) if presentation_mode else (20, 6))
    ax1 = fig1.add_subplot(axisbelow=True)
    ax2 = ax1.twinx() 
    ax1.plot(t, R / cpar.R_E, color = 'b', linewidth = linewidth)
    ax2.plot(t, T, color = 'r', linewidth = linewidth, linestyle = '-.')

    if num_sol.t[end_index] < 1e-3:
        ax1.set_xlabel('$t$ [μs]')
    else:
        ax1.set_xlabel('$t$ [ms]')
    ax1.set_ylabel('$R/R_E$ [-]', color = 'b')
    ax2.set_ylabel('$T$ [K]', color = 'r')
    if not presentation_mode: ax1.grid()
    
# textbox with initial conditions
    f"""Initial conditions:
    {'$R_E$':<25} {1e6*cpar.R_E: .2f}  $[\mu m]$
    {'$R_0/R_E$':<25} {cpar.ratio: .2f}  $[-]$
    {'$P_{amb}$':<25} {1e-5*cpar.P_amb: .2f}  $[bar]$
    {'$α_M$':<25} {cpar.alfa_M: .2f}  $[-]$
    {'$T_inf$':<25} {cpar.T_inf-273.15: .2f}  $[°C]$
    {'$P_{vapour}$':<25} {cpar.P_v: .1f}  $[Pa]$
    {'$μ_L$':<25} {1000*cpar.mu_L: .2f} $[mPa*s]$
    {'$surfactant$':<25} {cpar.surfactant: .2f}  $[-]$
    {'Initial content:':<20}
    """
    text = f'Initial conditions:\n'
    text += f'    $R_E$ = {1e6*cpar.R_E: .2f} $[\mu m]$\n'
    if cpar.ratio != 1.0:
        text += f'    $R_0/R_E$ = {cpar.ratio: .2f} $[-]$\n'
    text += f'    $P_{{amb}}$ = {1e-5*cpar.P_amb: .2f} $[bar]$\n'
    text += f'    $T_{{inf}}$ = {cpar.T_inf-273.15: .2f} $[°C]$\n'
    text += f'    $P_{{vapour}}$ = {cpar.P_v: .1f} $[Pa]$\n'
    text += f'Initial content:\n    '
    for gas, fraction in zip(cpar.gases, cpar.fractions):
        text += f'{int(100*fraction)}% {par.species[gas]}, ' 
    text = text[:-2] + f'\nExcitation = {excitation_type}:\n'
    for name, unit in zip(excitation_args, excitation_units):
        text += f'    {name} = {cpar[name]: .2f} [{unit}]\n'
    text = text[:-1]

    if show_cpar and not presentation_mode:
        ax1.text(
            0.98, 0.95, # coordinates
            text, transform=ax1.transAxes,
            horizontalalignment='right', verticalalignment='top', multialignment='left',
            fontsize=14, fontstyle='oblique',
            bbox={'facecolor': 'white', 'alpha': 1.0, 'pad': 10},
        )
    
    plt.show()

# plot reactions
    plt.rcParams.update({'font.size': 24 if presentation_mode else 18})
    fig2 = plt.figure(figsize=(16, 9) if presentation_mode else (20, 9))
    ax = fig2.add_subplot(axisbelow=True)

    # plot the lines
        # use this to generate colors:
            # import seaborn as sns
            # colors = sns.color_palette('Set1', n_colors=10)
            # print(colors.as_hex()); colors
    colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#d48044', '#33adff', '#a65628', '#f781bf', '#d444ca', '#d4ae44']
    color_index = 0
    texts = []
    max_mol = np.max(n, axis=1) # maximum amounts of species [mol]
    indexes_to_plot = np.argsort(max_mol)[-10:] if len(max_mol) >= 10 else np.argsort(max_mol) # Get the indexes of the 10 largest values
    for i, specie in enumerate(par.species):
        name = specie
        for digit in range(10): # turns 'H2O2' into 'H_2O_2'
            name = name.replace(str(digit), '_' + str(digit))
        if i in indexes_to_plot:
            color = colors[color_index]
            color_index = color_index + 1 if color_index < len(colors) - 1 else 0
            linewidth = 2.0 if presentation_mode and n[i, -1] > 1e-24 else 1.0
            ax.plot(t, n[i], linewidth = linewidth, color=color, label = '$' + name + '$') # PLOT HERE
            texts.append((color, name, n[i, -1]))            
        elif not presentation_mode:
            linewidth = 2.0 if presentation_mode else 1.0
            ax.plot(t, n[i], linewidth = linewidth, label = '$' + name + '$')  # PLOT HERE

    # make legend
    texts.sort(key=lambda x: x[-1], reverse=True)
    last_n_final = 1.0e100
    for text in texts:
        color, name, n_final = text
        # spaceing
        if n_final < 1e-24: continue
        limit = 5.5 if presentation_mode else 3.5
        if last_n_final / n_final < limit:
            n_final = last_n_final / limit
        last_n_final = n_final
        # place text
        ax.text(
            t[-1],
            n_final,
            '$' + name + '$',
            color=color,
            fontsize=24 if presentation_mode else 18,
            verticalalignment='center',
            bbox={'facecolor': 'white', 'pad': 0, 'linewidth': 0.0},
        )

    # plot settings
    plt.ylim([1e-24, 5.0*max_mol[indexes_to_plot[-1]]])
    ax.set_yscale('log')
    if num_sol.t[end_index] < 1e-3:
        ax.set_xlabel('$t$ [μs]')
    else:
        ax.set_xlabel('$t$ [ms]')
    ax.set_ylabel('$n_k$ [mol]')
    if not presentation_mode: ax.grid()
    if show_legend: ax.legend()

    plt.show()

# plot pressure excitation
    if plot_pressure:
        plt.rcParams.update({'font.size': 24 if presentation_mode else 18})
        fig3 = plt.figure(figsize=(16, 9) if presentation_mode else (20, 6))
        ax = fig3.add_subplot(axisbelow=True)
        ex_args = np.array([cpar[name] for name in excitation_args], dtype=np.float64)
        external_pressure = [1e-3*Excitation(time, cpar.P_amb, ex_args)[0] for time in num_sol.t[:end_index]] # [MPa]
        ax.plot(t, external_pressure, color='orange', label='external pressure', linewidth = 2.0 if presentation_mode else 1.0)

        if num_sol.t[end_index] < 1e-3:
            ax.set_xlabel('$t$ [μs]')
        else:
            ax.set_xlabel('$t$ [ms]')
        ax.set_ylabel('Pressure excitation [kPa]')
        if not presentation_mode: ax.grid()

        plt.show()

# plot pressure
    if plot_pressure:
        plt.rcParams.update({'font.size': 24 if presentation_mode else 18})
        fig4 = plt.figure(figsize=(16, 9) if presentation_mode else (20, 6))
        ax = fig4.add_subplot(axisbelow=True)
        ax.plot(t, internal_pressure, color = 'g', label='internal pressure', linewidth = 2.0 if presentation_mode else 1.0)

        if num_sol.t[end_index] < 1e-3:
            ax.set_xlabel('$t$ [μs]')
        else:
            ax.set_xlabel('$t$ [ms]')
        ax.set_ylabel('Internal pressure [MPa]')
        ax.set_yscale('log')
        if not presentation_mode: ax.grid()

        plt.show()
    
# saving the plots
    if base_name != '':
        try:
            metadata = {key: str(data[key]) for key in data.keys()}
            fig1.savefig(base_name + '_1.png', format='png', metadata=metadata)
            fig2.savefig(base_name + '_2.png', format='png', metadata=metadata)
            if plot_pressure:
                fig3.savefig(base_name + '_3.png', format='png', metadata=metadata)
                fig4.savefig(base_name + '_4.png', format='png', metadata=metadata)
        except:
            print(print(colored(f'Error in saving {base_name}_1.png','red')))

# print data
    print_data(cpar, data)
    return None
           

"""________________________________Save to CSV________________________________"""

def get_settings_and_info():
    return f'''General information:
    Created: {datetime.now().strftime("%Y.%m.%d %H:%M:%S (YYYY.MM.DD HH:MM:SS)")}
    Computer name: {socket.gethostname()}
    User name: {os.getlogin()}
    Number of cores: {psutil.cpu_count(logical=False)}
    Number of logical threads: {psutil.cpu_count(logical=True)}
    CPU frequency: {psutil.cpu_freq().max: .2f} MHz
    RAM size: {psutil.virtual_memory().total / 1024**2: .2f} MB

parameters settings:
    model = {par.model}
    species = {par.species}
    number of species = {par.K}
    number of reactions = {par.I}

full_bubble_model settings:
    enable_heat_transfer = {enable_heat_transfer}
    enable_evaporation = {enable_evaporation} 
    enable_reactions = {enable_reactions}
    enable_dissipated_energy = {enable_dissipated_energy}
    target_specie = \'{target_specie}\' # Specie to calculate energy effiqiency
    excitation_type = \'{excitation_type}\' # function to calculate pressure excitation
'''

class Make_dir:
    # constructor
    def __init__(self, folder_name, file_base_name='output_', separator=','):
        self.folder_name = folder_name
        self.file_base_name = file_base_name
        self.separator = separator
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
            element = str(element).replace(',', ' ').replace('[', '').replace(']', '')
            if isinstance(element, float):
                line += f'{float(element): e}' + self.separator
            else:     
                line += element + self.separator
        return line[:-1]
        
    # writes a data dict into the currently opened file
    def write_line(self, data):
        line = self.list_to_string([data[key] for key in keys])
        line += self.separator + self.list_to_string([x for x in data['x_initial'][:-1]] + [x for x in data['x_final'][:-1]])
        self.file.write(line + '\n')
        self.lines += 1
        
    # saves a numerical solution
    def write_solution(self, data, num_sol, file_base_name):
        # create file containing data
        file = os.path.join(self.save_dir, file_base_name + '_data.csv')
        file = open(file, 'w')
        # write header line
        line = self.list_to_string(keys + ['R_0', 'R_dot_0', 'T_0'] + ['c_' + specie + '_0' for specie in par.species] + ['R_last', 'R_dot_last', 'T_last'] + ['c_' + specie + '_last' for specie in par.species])
        file.write(line + '\n')
        # write data
        line = self.list_to_string([data[key]] for key in keys)
        line += self.separator + self.list_to_string([x for x in data.x_initial[:-1]] + [x for x in data.x_final[:-1]] + [data.expansion_work,data.dissipated_acoustic_energy,data.energy_efficiency])
        file.write(line + '\n')
        file.close()

    # create file containing num_sol
        file = os.path.join(self.save_dir, file_base_name + '_num_sol.csv')
        file = open(file, 'w')
        # write header line
        line = self.list_to_string(['t', 'R', 'R_dot', 'T'] + ['c_' + specie for specie in par.species] + ['dissipated_acoustic_energy'])
        file.write(line + '\n')
        # write data
        for i in range(len(num_sol.t)):
            line = self.list_to_string([num_sol.t[i]] + list(num_sol.y[:, i]))
            file.write(line + '\n')
        file.close()
    
    # save any string
    def write_string(self, string, file_base_name):
        # create file
        file = os.path.join(self.save_dir, file_base_name + '.txt')
        try:
            file = open(file, 'x')
        except FileExistsError:
            print(colored(f'Error, file \'{file}\' already exists. ', 'red'))
            return None
        # create header
        file.write(get_settings_and_info())
        # write string
        file.write(string)
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
        line = self.list_to_string(keys + ['R_0', 'R_dot_0', 'T_0'] + ['c_' + specie + '_0' for specie in par.species] + ['R_last', 'R_dot_last', 'T_last'] + ['c_' + specie + '_last' for specie in par.species])
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