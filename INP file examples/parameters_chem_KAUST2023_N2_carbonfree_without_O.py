"""
This is an automatically generated file.
This python file contains all the numerical constants and coefficents required for the simulations.
Recommended usage:
   importing: import parameters as par
   usage: print('Universal gas constant =', par.R_g)
"""

model = 'chem_KAUST2023_N2_carbonfree_without_O'
import numpy as np

"""________________________________Physical constants________________________________"""

c_L            = 1483.0                   # Liquid sound speed at 30 °C [m/s]
rho_L          = 998.2                    # Liquid density [kg/m^3]
C_p_L          = 4178.0                   # Isobar heat capacity of water [J/(kg*K)]
sigma          = 0.07197                  # Surface tension [N/m]
mu_L           = 0.001                    # Dynamic viscosity at 30 °C and 1 atm [Pa*s]
P_v            = 2338.1                   # Saturated vapour pressure at 30 °C [Pa]
alfa_M         = 0.35                     # water accommodation coefficient [-]
R_g            = 8.31446                  # Universal gas constant [J/mol/K]
R_erg          = 83144600.0               # Universal gas constant [erg/mol/K]
R_cal          = 1.987204                 # Universal gas constant [cal/mol/K]
N_A            = 6.02214e+23              # Avogadro's number [-]
h              = 6.62607015e-34           # Planck constant [m^2*kg/s]
k_B            = 1.3806487394846352e-23   # Boltzmann constant [J/K]
R_v            = 461.521126               # Specific gas constant of water [J/kg/K]
erg2J          = 1e-07                    # Conversion factor from erg to J
J2erg          = 10000000.0               # Conversion factor from J to erg
cal2J          = 4.184                    # Conversion factor from cal to J
J2cal          = 0.2390057361376673       # Conversion factor from J to cal
atm2Pa         = 101325.0                 # Conversion factor from atm to Pa
bar2Pa         = 100000.0                 # Conversion factor from bar to Pa
absolute_zero  = 273.15                   # Zero °C in Kelvin


"""________________________________Species________________________________"""

elements = np.array(['AR','C','H','HE','N'])
#                            0,         1,         2,         3,         4,         5,         6,         7,         8
species = np.array([      'AR',       'H',      'H2',    'H2NN',       'N',      'N2',    'N2H2',    'N2H3',    'N2H4'])

# molar mass [g/mol]
W = np.array([          39.948,   1.00797,   2.01594,  30.02934,   14.0067,   28.0134,  30.02934,  31.03731,  32.04528], dtype=np.float64)

# thermal conductivity [W / m / K]
lambdas = np.array([    0.0177,       0.0,    0.1805,       0.0,       0.0,   0.02583,       0.0,       0.0,       0.0], dtype=np.float64)

index = dict(
        AR= 0,      H= 1,     H2= 2,   H2NN= 3,      N= 4,     N2= 5,   N2H2= 6,   N2H3= 7,   N2H4= 8
)

indexOfWater = -1
K = 9   # Number of species


"""________________________________NASA polynomials________________________________"""

N = 5    # degree of polynomials
TempRange = np.array([
    #   T_low   T_high    T_mid 
    [   300.0,  5000.0,  1000.0],    # AR
    [   200.0,  3500.0,  1000.0],    # H
    [   200.0,  3500.0,  1000.0],    # H2
    [   200.0,  6000.0,  1000.0],    # H2NN
    [   300.0,  5000.0,  1000.0],    # N
    [   300.0,  5000.0,  1000.0],    # N2
    [  298.15,  2000.0,  1000.0],    # N2H2
    [   200.0,  6000.0,  1000.0],    # N2H3
    [   200.0,  6000.0,  1000.0]     # N2H4
], dtype=np.float64)

# LOW NASA coefficients
a_low = np.array([
    #             a_1              a_2              a_3              a_4              a_5              a_6              a_7 
    [             2.5,             0.0,             0.0,             0.0,             0.0,        -745.375,           4.366],    # AR
    [             2.5,  7.05332819e-13, -1.99591964e-15,  2.30081632e-18, -9.27732332e-22,      25473.6599,    -0.446682853],    # H
    [      2.34433112,   0.00798052075,  -1.9478151e-05,  2.01572094e-08, -7.37611761e-12,     -917.935173,     0.683010238],    # H2
    [      4.53204001,  -0.00732418578,  3.00803713e-05, -3.04000551e-08,  1.04700639e-11,      34958.0003,      1.51074195],    # H2NN
    [        2.503071,   -2.180018e-05,    5.420529e-08,    -5.64756e-11,    2.099904e-14,         56098.9,        4.167566],    # N
    [        3.298677,    0.0014082404,   -3.963222e-06,    5.641515e-09,   -2.444854e-12,      -1020.8999,        3.950372],    # N2
    [      4.14915878,  -0.00482013543,  2.19606525e-05, -2.03662902e-08,    6.253672e-12,      23287.4597,      3.19154222],    # N2H2
    [      3.42125505,    0.0013490159,  2.23459071e-05, -2.99727732e-08,   1.2097897e-11,      25305.6139,      7.83176309],    # N2H3
    [      3.83472149, -0.000649129555,  3.76848463e-05, -5.00709182e-08,  2.03362064e-11,      10089.3925,       5.7527203]     # N2H4
], dtype=np.float64)

# LOW NASA coefficients
a_high = np.array([
    #             a_1              a_2              a_3              a_4              a_5              a_6              a_7 
    [             2.5,             0.0,             0.0,             0.0,             0.0,        -745.375,           4.366],    # AR
    [      2.50000001, -2.30842973e-11,  1.61561948e-14, -4.73515235e-18,  4.98197357e-22,      25473.6599,    -0.446682914],    # H
    [       3.3372792, -4.94024731e-05,  4.99456778e-07, -1.79566394e-10,  2.00255376e-14,     -950.158922,     -3.20502331],    # H2
    [       3.0590367,   0.00618382347, -2.22171165e-06,  3.58539206e-10, -2.14532905e-14,      34853.0149,      6.69893515],    # H2NN
    [        2.450268,    0.0001066146,   -7.465337e-08,    1.879652e-11,   -1.025984e-15,        56116.04,        4.448758],    # N
    [         2.92664,    0.0014879768,    -5.68476e-07,   1.0097038e-10,   -6.753351e-15,       -922.7977,        5.980528],    # N2
    [       2.3838838,   0.00663619559, -1.81669068e-06, -1.19497031e-10,  9.31659255e-14,      23420.7532,      10.6091658],    # N2H2
    [      4.04483566,   0.00731130186, -2.47625799e-06,  3.83733021e-10, -2.23107573e-14,      24809.8603,      2.88423392],    # N2H3
    [      4.93957357,   0.00875017187, -2.99399058e-06,  4.67278418e-10, -2.73068599e-14,      9282.65548,     -2.69439772]     # N2H4
], dtype=np.float64)


"""________________________________Reaction constants________________________________"""

I = 49    # Number of reactions
# Pre-exponential factors [cm^3/mol/s v 1/s]
A = np.array([
               4.577e+19,            5.84e+18,            5.84e+18,            1.89e+18,      165000000000.0,
         6260000000000.0,    56340000000000.0,    30000000000000.0,  1200000000000000.0,            109000.0,
        69000000000000.0,                0.57,               5.636,              9574.0,             1.9e+27,
             482000000.0,           2400000.0,            271000.0,   130000000000000.0,         170000000.0,
       560000000000000.0,            5.19e+38,            276000.0,                0.76,       25000000000.0,
         1000000000000.0,      128000000000.0,              7476.0,           6243000.0,      100000000000.0,
                   0.608,                11.1,    20000000000000.0,    10000000000000.0,     3000000000000.0,
        10000000000000.0,             72000.0,             3.4e+26,         480000000.0,    70000000000000.0,
               1800000.0,             2.2e+16,           2890000.0,      100000000000.0,         330000000.0,
       100000000000000.0,    50000000000000.0,    50000000000000.0,    10000000000000.0
], dtype=np.float64)

# Temperature exponent [-]
b = np.array([
                    -1.4,                -1.1,                -1.1,               -0.85,                0.71,
                  -0.036,              -0.036,                 0.0,                 0.0,                2.59,
                     0.0,                3.88,                3.53,                2.46,               -3.05,
                    1.76,                 2.0,               2.226,              -0.272,                1.62,
                  -0.414,               -7.84,                2.56,                 4.0,                 0.5,
                     0.0,               0.819,               2.796,                1.89,                 0.0,
                   3.574,                3.08,                 0.0,                 0.0,                 0.0,
                     0.0,                1.88,               -4.83,                 1.5,                 0.0,
                    1.94,                 0.0,                2.23,                 0.5,                 0.0,
                     0.0,                 0.0,                 0.0,                 0.0
], dtype=np.float64)

# Activation energy [cal/mol]
E = np.array([
            104379.78643,        104379.78643,        104379.78643,        224999.52768,            931.0032,
              -160.90358,          -160.90358,                 0.0,          76002.5903,           1812.3264,
                     0.0,           341.99712,           552.60058,            107.3088,         66099.85603,
                739.2384,         -1192.00205,          -1030.0055,             -77.004,         11782.98317,
                65.99491,         67099.85482,           1220.0017,          4047.98602,         29799.93197,
              1990.00195,         48099.89779,           4679.9951,           247.00896,                 0.0,
              1189.99498,            211.0009,                 0.0,           9939.9744,                 0.0,
              9934.98653,          8801.98445,         46219.90723,          -894.00154,                 0.0,
             -1151.99971,          93469.8001,         10399.97146,         21599.94989,                 0.0,
                     0.0,                 0.0,                 0.0,          9934.98653
], dtype=np.float64)


"""________________________________Reaction matrixes________________________________"""

# Forward reaction matrix
nu_forward = np.array([
    #  AR    H   H2 H2NN    N   N2 N2H2 N2H3 N2H4 
    [   0,   0,   1,   0,   0,   0,   0,   0,   0],    #  0. H2+M=2H+M
    [   1,   0,   1,   0,   0,   0,   0,   0,   0],    #  1. H2+AR=2H+AR
    [   0,   0,   1,   0,   0,   0,   0,   0,   0],    #  2. H2+HE=2H+HE
    [   0,   0,   0,   0,   0,   1,   0,   0,   0],    #  3. N2+M=2N+M
    [   0,   1,   0,   0,   0,   0,   0,   0,   0],    #  4. NH+H=N+H2
    [   0,   0,   0,   0,   0,   0,   0,   0,   0],    #  5. 2NH=>N2+H2
    [   0,   0,   0,   0,   0,   0,   0,   0,   0],    #  6. 2NH=>N2+2H
    [   0,   0,   0,   0,   1,   0,   0,   0,   0],    #  7. NH+N=N2+H
    [   0,   0,   0,   0,   0,   0,   0,   0,   0],    #  8. NH2+M=NH+H+M
    [   0,   1,   0,   0,   0,   0,   0,   0,   0],    #  9. NH2+H=NH+H2
    [   0,   0,   0,   0,   1,   0,   0,   0,   0],    # 10. NH2+N=N2+2H
    [   0,   0,   0,   0,   0,   0,   0,   0,   0],    # 11. 2NH=NH2+N
    [   0,   0,   0,   0,   0,   0,   0,   0,   0],    # 12. 2NH2=NH3+NH
    [   0,   0,   0,   0,   0,   0,   0,   0,   0],    # 13. NH2+NH=NH3+N
    [   0,   0,   0,   0,   0,   0,   1,   0,   0],    # 14. N2H2+M=NNH+H+M
    [   0,   1,   0,   0,   0,   0,   1,   0,   0],    # 15. N2H2+H=NNH+H2
    [   0,   0,   0,   0,   0,   0,   1,   0,   0],    # 16. N2H2+NH=NNH+NH2
    [   0,   0,   0,   0,   0,   0,   1,   0,   0],    # 17. N2H2+NH2=NNH+NH3
    [   0,   0,   0,   0,   0,   0,   0,   0,   0],    # 18. NH2+NH=N2H2+H
    [   0,   0,   0,   0,   0,   0,   0,   0,   0],    # 19. 2NH2=N2H2+H2
    [   0,   0,   0,   0,   0,   0,   0,   0,   0],    # 20. 2NH2(+M)=N2H4(+M)
    [   0,   0,   0,   0,   0,   0,   0,   0,   1],    # 21. N2H4=N2H2+H2
    [   0,   1,   0,   0,   0,   0,   0,   0,   1],    # 22. N2H4+H=N2H3+H2
    [   0,   0,   0,   0,   0,   0,   0,   0,   1],    # 23. N2H4+NH2=N2H3+NH3
    [   0,   0,   0,   0,   0,   0,   1,   0,   1],    # 24. N2H4+N2H2=2N2H3
    [   0,   0,   0,   0,   0,   0,   0,   0,   1],    # 25. N2H4+NH=N2H3+NH2
    [   0,   0,   0,   0,   0,   0,   0,   1,   0],    # 26. N2H3(+M)=N2H2+H(+M)
    [   0,   1,   0,   0,   0,   0,   0,   1,   0],    # 27. N2H3+H=N2H2+H2
    [   0,   1,   0,   0,   0,   0,   0,   1,   0],    # 28. N2H3+H=H2NN+H2
    [   0,   1,   0,   0,   0,   0,   0,   1,   0],    # 29. N2H3+H=NH3+NH
    [   0,   0,   0,   0,   0,   0,   0,   1,   0],    # 30. N2H3+NH2=N2H2+NH3
    [   0,   0,   0,   0,   0,   0,   0,   1,   0],    # 31. N2H3+NH2=H2NN+NH3
    [   0,   0,   0,   0,   0,   0,   0,   1,   0],    # 32. N2H3+NH=N2H2+NH2
    [   0,   0,   0,   0,   0,   0,   1,   1,   0],    # 33. N2H3+N2H2=N2H4+NNH
    [   0,   0,   0,   0,   0,   0,   0,   2,   0],    # 34. 2N2H3=2NH3+N2
    [   0,   0,   0,   0,   0,   0,   2,   0,   0],    # 35. 2N2H2=N2H3+NNH
    [   0,   0,   0,   0,   0,   0,   0,   0,   0],    # 36. 2NH2=H2NN+H2
    [   0,   0,   0,   1,   0,   0,   0,   0,   0],    # 37. H2NN=NNH+H
    [   0,   1,   0,   1,   0,   0,   0,   0,   0],    # 38. H2NN+H=NNH+H2
    [   0,   1,   0,   1,   0,   0,   0,   0,   0],    # 39. H2NN+H=N2H2+H
    [   0,   0,   0,   1,   0,   0,   0,   0,   0],    # 40. H2NN+NH2=NNH+NH3
    [   0,   0,   0,   0,   0,   0,   0,   0,   0],    # 41. NH3+M=NH2+H+M
    [   0,   1,   0,   0,   0,   0,   0,   0,   0],    # 42. NH3+H=NH2+H2
    [   0,   0,   0,   0,   0,   0,   0,   0,   0],    # 43. NH3+NH2=N2H3+H2
    [   0,   0,   0,   0,   0,   0,   0,   0,   0],    # 44. NNH=N2+H
    [   0,   1,   0,   0,   0,   0,   0,   0,   0],    # 45. NNH+H=N2+H2
    [   0,   0,   0,   0,   0,   0,   0,   0,   0],    # 46. NNH+NH=N2+NH2
    [   0,   0,   0,   0,   0,   0,   0,   0,   0],    # 47. NNH+NH2=N2+NH3
    [   0,   0,   0,   0,   0,   0,   0,   0,   0]     # 48. 2NNH=N2H2+N2
], dtype=np.float64)

# Backward reaction matrix
nu_backward = np.array([
    #  AR    H   H2 H2NN    N   N2 N2H2 N2H3 N2H4 
    [   0,   2,   0,   0,   0,   0,   0,   0,   0],    #  0. H2+M=2H+M
    [   1,   2,   0,   0,   0,   0,   0,   0,   0],    #  1. H2+AR=2H+AR
    [   0,   2,   0,   0,   0,   0,   0,   0,   0],    #  2. H2+HE=2H+HE
    [   0,   0,   0,   0,   2,   0,   0,   0,   0],    #  3. N2+M=2N+M
    [   0,   0,   1,   0,   1,   0,   0,   0,   0],    #  4. NH+H=N+H2
    [   0,   0,   1,   0,   0,   1,   0,   0,   0],    #  5. 2NH=>N2+H2
    [   0,   2,   0,   0,   0,   1,   0,   0,   0],    #  6. 2NH=>N2+2H
    [   0,   1,   0,   0,   0,   1,   0,   0,   0],    #  7. NH+N=N2+H
    [   0,   1,   0,   0,   0,   0,   0,   0,   0],    #  8. NH2+M=NH+H+M
    [   0,   0,   1,   0,   0,   0,   0,   0,   0],    #  9. NH2+H=NH+H2
    [   0,   2,   0,   0,   0,   1,   0,   0,   0],    # 10. NH2+N=N2+2H
    [   0,   0,   0,   0,   1,   0,   0,   0,   0],    # 11. 2NH=NH2+N
    [   0,   0,   0,   0,   0,   0,   0,   0,   0],    # 12. 2NH2=NH3+NH
    [   0,   0,   0,   0,   1,   0,   0,   0,   0],    # 13. NH2+NH=NH3+N
    [   0,   1,   0,   0,   0,   0,   0,   0,   0],    # 14. N2H2+M=NNH+H+M
    [   0,   0,   1,   0,   0,   0,   0,   0,   0],    # 15. N2H2+H=NNH+H2
    [   0,   0,   0,   0,   0,   0,   0,   0,   0],    # 16. N2H2+NH=NNH+NH2
    [   0,   0,   0,   0,   0,   0,   0,   0,   0],    # 17. N2H2+NH2=NNH+NH3
    [   0,   1,   0,   0,   0,   0,   1,   0,   0],    # 18. NH2+NH=N2H2+H
    [   0,   0,   1,   0,   0,   0,   1,   0,   0],    # 19. 2NH2=N2H2+H2
    [   0,   0,   0,   0,   0,   0,   0,   0,   1],    # 20. 2NH2(+M)=N2H4(+M)
    [   0,   0,   1,   0,   0,   0,   1,   0,   0],    # 21. N2H4=N2H2+H2
    [   0,   0,   1,   0,   0,   0,   0,   1,   0],    # 22. N2H4+H=N2H3+H2
    [   0,   0,   0,   0,   0,   0,   0,   1,   0],    # 23. N2H4+NH2=N2H3+NH3
    [   0,   0,   0,   0,   0,   0,   0,   2,   0],    # 24. N2H4+N2H2=2N2H3
    [   0,   0,   0,   0,   0,   0,   0,   1,   0],    # 25. N2H4+NH=N2H3+NH2
    [   0,   1,   0,   0,   0,   0,   1,   0,   0],    # 26. N2H3(+M)=N2H2+H(+M)
    [   0,   0,   1,   0,   0,   0,   1,   0,   0],    # 27. N2H3+H=N2H2+H2
    [   0,   0,   1,   1,   0,   0,   0,   0,   0],    # 28. N2H3+H=H2NN+H2
    [   0,   0,   0,   0,   0,   0,   0,   0,   0],    # 29. N2H3+H=NH3+NH
    [   0,   0,   0,   0,   0,   0,   1,   0,   0],    # 30. N2H3+NH2=N2H2+NH3
    [   0,   0,   0,   1,   0,   0,   0,   0,   0],    # 31. N2H3+NH2=H2NN+NH3
    [   0,   0,   0,   0,   0,   0,   1,   0,   0],    # 32. N2H3+NH=N2H2+NH2
    [   0,   0,   0,   0,   0,   0,   0,   0,   1],    # 33. N2H3+N2H2=N2H4+NNH
    [   0,   0,   0,   0,   0,   1,   0,   0,   0],    # 34. 2N2H3=2NH3+N2
    [   0,   0,   0,   0,   0,   0,   0,   1,   0],    # 35. 2N2H2=N2H3+NNH
    [   0,   0,   1,   1,   0,   0,   0,   0,   0],    # 36. 2NH2=H2NN+H2
    [   0,   1,   0,   0,   0,   0,   0,   0,   0],    # 37. H2NN=NNH+H
    [   0,   0,   1,   0,   0,   0,   0,   0,   0],    # 38. H2NN+H=NNH+H2
    [   0,   1,   0,   0,   0,   0,   1,   0,   0],    # 39. H2NN+H=N2H2+H
    [   0,   0,   0,   0,   0,   0,   0,   0,   0],    # 40. H2NN+NH2=NNH+NH3
    [   0,   1,   0,   0,   0,   0,   0,   0,   0],    # 41. NH3+M=NH2+H+M
    [   0,   0,   1,   0,   0,   0,   0,   0,   0],    # 42. NH3+H=NH2+H2
    [   0,   0,   1,   0,   0,   0,   0,   1,   0],    # 43. NH3+NH2=N2H3+H2
    [   0,   1,   0,   0,   0,   1,   0,   0,   0],    # 44. NNH=N2+H
    [   0,   0,   1,   0,   0,   1,   0,   0,   0],    # 45. NNH+H=N2+H2
    [   0,   0,   0,   0,   0,   1,   0,   0,   0],    # 46. NNH+NH=N2+NH2
    [   0,   0,   0,   0,   0,   1,   0,   0,   0],    # 47. NNH+NH2=N2+NH3
    [   0,   0,   0,   0,   0,   1,   1,   0,   0]     # 48. 2NNH=N2H2+N2
], dtype=np.float64)

nu = nu_backward - nu_forward

reaction_order = np.array([
    [   1],    #  0. H2+M=2H+M
    [   2],    #  1. H2+AR=2H+AR
    [   1],    #  2. H2+HE=2H+HE
    [   1],    #  3. N2+M=2N+M
    [   1],    #  4. NH+H=N+H2
    [   0],    #  5. 2NH=>N2+H2
    [   0],    #  6. 2NH=>N2+2H
    [   1],    #  7. NH+N=N2+H
    [   0],    #  8. NH2+M=NH+H+M
    [   1],    #  9. NH2+H=NH+H2
    [   1],    # 10. NH2+N=N2+2H
    [   0],    # 11. 2NH=NH2+N
    [   0],    # 12. 2NH2=NH3+NH
    [   0],    # 13. NH2+NH=NH3+N
    [   1],    # 14. N2H2+M=NNH+H+M
    [   2],    # 15. N2H2+H=NNH+H2
    [   1],    # 16. N2H2+NH=NNH+NH2
    [   1],    # 17. N2H2+NH2=NNH+NH3
    [   0],    # 18. NH2+NH=N2H2+H
    [   0],    # 19. 2NH2=N2H2+H2
    [   0],    # 20. 2NH2(+M)=N2H4(+M)
    [   1],    # 21. N2H4=N2H2+H2
    [   2],    # 22. N2H4+H=N2H3+H2
    [   1],    # 23. N2H4+NH2=N2H3+NH3
    [   2],    # 24. N2H4+N2H2=2N2H3
    [   1],    # 25. N2H4+NH=N2H3+NH2
    [   1],    # 26. N2H3(+M)=N2H2+H(+M)
    [   2],    # 27. N2H3+H=N2H2+H2
    [   2],    # 28. N2H3+H=H2NN+H2
    [   2],    # 29. N2H3+H=NH3+NH
    [   1],    # 30. N2H3+NH2=N2H2+NH3
    [   1],    # 31. N2H3+NH2=H2NN+NH3
    [   1],    # 32. N2H3+NH=N2H2+NH2
    [   2],    # 33. N2H3+N2H2=N2H4+NNH
    [   2],    # 34. 2N2H3=2NH3+N2
    [   2],    # 35. 2N2H2=N2H3+NNH
    [   0],    # 36. 2NH2=H2NN+H2
    [   1],    # 37. H2NN=NNH+H
    [   2],    # 38. H2NN+H=NNH+H2
    [   2],    # 39. H2NN+H=N2H2+H
    [   1],    # 40. H2NN+NH2=NNH+NH3
    [   0],    # 41. NH3+M=NH2+H+M
    [   1],    # 42. NH3+H=NH2+H2
    [   0],    # 43. NH3+NH2=N2H3+H2
    [   0],    # 44. NNH=N2+H
    [   1],    # 45. NNH+H=N2+H2
    [   0],    # 46. NNH+NH=N2+NH2
    [   0],    # 47. NNH+NH2=N2+NH3
    [   0]     # 48. 2NNH=N2H2+N2
], dtype=np.int64)


"""________________________________Three-body reactions________________________________"""

ThirdBodyIndexes = np.array([   0,   3,   8,  14,  20,  26,  41], dtype=np.int64)
ThirdBodyCount = 7

# third-body efficiency factors
alfa = np.array([
    #      AR        H       H2     H2NN        N       N2     N2H2     N2H3     N2H4 
    [     0.0,     1.0,     2.5,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0],    #  0. H2+M=2H+M
    [     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0],    #  3. N2+M=2N+M
    [     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0],    #  8. NH2+M=NH+H+M
    [     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0],    # 14. N2H2+M=NNH+H+M
    [    0.59,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0],    # 20. 2NH2(+M)=N2H4(+M)
    [     1.0,     1.0,     1.0,     1.0,     1.0,     2.0,     1.0,     1.0,     1.0],    # 26. N2H3(+M)=N2H2+H(+M)
    [     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0]     # 41. NH3+M=NH2+H+M
], dtype=np.float64)


"""________________________________Irreversible reactions________________________________"""

IrreversibleIndexes = np.array([   5,   6], dtype=np.int64)
IrreversibleCount = 2


"""________________________________Pressure-dependent reactions________________________________"""

PressureDependentIndexes = np.array([  20,  26], dtype=np.int64)
PressureDependentCount = 2

LindemannIndexes = np.array([], dtype=np.int64)
LindemannCount = 0

# Fall-off parameters
ReacConst = np.array([
    #               A_0                b_0                E_0 
    [           1.6e+34,             -5.49,        1987.00128],    # 20. 2NH2(+M)=N2H4(+M)
    [          3.84e+40,             -6.88,      54499.893984]     # 26. N2H3(+M)=N2H2+H(+M)
], dtype=np.float64)

TroeIndexes = np.array([  20,  26], dtype=np.int64)
TroeCount = 2

# Troe parameters
Troe = np.array([
    #              alfa               T***                 T*                T** 
    [              0.31,             1e-30,             1e+30,             1e+30],    # 20. 2NH2(+M)=N2H4(+M)
    [             0.842,           80000.0,              28.0,            7300.0]     # 26. N2H3(+M)=N2H2+H(+M)
], dtype=np.float64)

SRIIndexes = np.array([], dtype=np.int64)
SRICount = 0

# SRI parameters
SRI = np.array([
    #                 a                  b                  c                  d                  e 
    [] 
], dtype=np.float64)

PlogIndexes = np.array([], dtype=np.int64)
PlogStart = np.array([   0], dtype=np.int64)
PlogStop = np.array([   1], dtype=np.int64)
PlogCount = 0

# PLOG parameters
Plog = np.array([
    #               P_1                A_1                b_1                E_1 
    [] 
], dtype=np.float64)

