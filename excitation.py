"""
This file contains different pressure excitation functions.
Specify which one to use under "Settings" in full_bubble_model.py
"""

import numpy as np
from termcolor import colored
from numba import njit   # Just In Time compiler
from numba.types import float64   # JIT types

def getExcitation(excitation_type='no_excitation'):
    """
    Returns the pressure excitation function and the arguments it takes.
    Available excitation types:
     * 'no_excitation': constant ambient pressure
     * 'two_sinusoids': two sinusoids with different frequencies and amplitudes, and a phase shift between them
     * 'sin_sqr': sinusoid squared with only n amplitude cycles
     * 'slow_expansion': decrease from abient pressure to min_pressure (decay_time), followed by a slow increase back to ambient pressure (increase_time)
     * 'sin_impulse_flat_ends': sinusoid with only n amplitude cycles, the ends are smoothed out
     * 'sin_impulse': sinusoid with only n amplitude cycles, the ends are not smoothed out
    """

    if excitation_type == 'no_excitation':
        @njit(float64[:](float64, float64, float64[:]))
        def Excitation(t, P_amb, args):
            P_amb = args[0]
            return np.array([P_amb, 0.0], dtype=float64)
        
        args = []
        units = []
        return Excitation, args, units
    
    elif excitation_type == 'two_sinusoids':
        @njit(float64[:](float64, float64, float64[:]))
        def Excitation(t, P_amb, args):
            p_A1, p_A2, freq1, freq2, theta_phase = args
            p_Inf = P_amb + p_A1*np.sin(2.0*np.pi*freq1*t) + p_A2*np.sin(2.0*np.pi*freq2*t + theta_phase) 
            p_Inf_dot = p_A1*2.0*np.pi*freq1*np.cos(2.0*np.pi*freq1*t) + p_A2*2.0*np.pi*freq2*np.cos(2.0*np.pi*freq2*t+theta_phase)
            return np.array([p_Inf, p_Inf_dot], dtype=float64)
        
        args = ['p_A1', 'p_A2', 'freq1', 'freq2', 'theta_phase']
        units = ['Pa', 'Pa', 'Hz', 'Hz', 'rad']
        return Excitation, args, units
    
    elif excitation_type == 'sin_sqr':
        @njit(float64[:](float64, float64, float64[:]))
        def Excitation(t, P_amb, args):
            p_A, freq, n = args
            if t < n / freq:
                p_Inf = P_amb + p_A * np.sin(2.0*np.pi*freq*t)**2
                p_Inf_dot = 4.0*p_A*np.pi*freq * np.sin(4.0*np.pi*freq*t)
            else:
                p_Inf = P_amb
                p_Inf_dot = 0.0
            return np.array([p_Inf, p_Inf_dot], dtype=np.float64)
        
        args = ['p_A', 'freq', 'n']
        units = ['Pa', 'Hz', '-']
        return Excitation, args, units
    
    elif excitation_type == 'slow_expansion':
        @njit(float64[:](float64, float64, float64[:]))
        def Excitation(t, P_amb, args):
            decay_time, increase_time, min_pressure = args
            if t < 0.0:
                p_Inf = P_amb
                p_Inf_dot = 0.0
            elif t > decay_time+increase_time:
                p_Inf = P_amb
                p_Inf_dot = 0.0
            elif t < decay_time:
                p_Inf = P_amb - (P_amb-min_pressure) / decay_time * t
                p_Inf_dot = - (P_amb-min_pressure) / decay_time
            else:
                p_Inf = min_pressure + (P_amb-min_pressure) / (increase_time) * (t-decay_time)
                p_Inf_dot = (P_amb-min_pressure) / (increase_time)

            return np.array([p_Inf, p_Inf_dot], dtype=np.float64)
        
        args = ['decay_time', 'increase_time', 'min_pressure']
        units = ['s', 's', 'Hz']
        return Excitation, args, units
    
    elif excitation_type == 'sin_impulse_flat_ends':
        @njit(float64[:](float64, float64, float64[:]))
        def Excitation(t, P_amb, args):
            p_A, freq, n = args
            if t < 0.0:
                p_Inf = P_amb
                p_Inf_dot = 0.0
            elif t > n / freq:
                p_Inf = P_amb
                p_Inf_dot = 0.0
            elif t < 0.25 / freq:
                p_Inf = P_amb + p_A*np.sin(2.0*np.pi*freq*t)**2
                p_Inf_dot = 4.0*p_A*np.pi*freq * np.sin(4.0*np.pi*freq*t)
            elif t > (n - 0.25) / freq:
                sign = np.sign( np.sin(2.0*np.pi*freq*t) )
                p_Inf = P_amb + p_A*np.sin(2.0*np.pi*freq*t)**2 * sign
                p_Inf_dot = 4.0*p_A*np.pi*freq * np.sin(4.0*np.pi*freq*t) * sign
            else:
                p_Inf = P_amb + p_A*np.sin(2.0*np.pi*freq*t)
                p_Inf_dot = p_A*2.0*np.pi*freq*np.cos(2.0*np.pi*freq*t)

            return np.array([p_Inf, p_Inf_dot], dtype=np.float64)
        
        args = ['p_A', 'freq', 'n']
        units = ['Pa', 'Hz', '-']
        return Excitation, args, units
    
    elif excitation_type == 'sin_impulse':
        @njit(float64[:](float64, float64, float64[:]))
        def Excitation(t, P_amb, args):
            p_A, freq, n = args
            if t < 0.0:
                p_Inf = P_amb
                p_Inf_dot = 0.0
            elif t > n / freq:
                p_Inf = P_amb
                p_Inf_dot = 0.0
            else:
                insin = 2.0*np.pi*freq
                p_Inf = P_amb + p_A*np.sin(insin*t)
                p_Inf_dot = p_A*insin*np.cos(insin*t)

            return np.array([p_Inf, p_Inf_dot], dtype=np.float64)
        
        args = ['p_A', 'freq', 'n']
        units = ['Pa', 'Hz', '-']
        return Excitation, args, units
    
    else:
        print(colored(f'Warning: Excitation excitation_type \'{excitation_type}\' not recognized. Using \'no_excitation\' instead. ', 'yellow'))
        return getExcitation('no_excitation')
    