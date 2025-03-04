"""
This is an automatically generated file.
This python file contains all the numerical constants and coefficents required for the simulations.
Recommended usage:
   importing: import parameters as par
   usage: print('Universal gas constant =', par.R_g)
"""

model = 'chemkin_AR_HE_FIXED_by_Cantera'
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

elements = np.array(['C','H','N','O','AR','HE'])
#                            0,         1,         2,         3,         4,         5,         6,         7,         8,         9,        10,        11
species = np.array([       'H',      'H2',       'O',      'O2',      'OH',     'H2O',      'N2',     'HO2',    'H2O2',      'AR',      'HE',    'OHEX'])

# molar mass [g/mol]
W = np.array([         1.00797,   2.01594,   15.9994,   31.9988,  17.00737,  18.01534,   28.0134,  33.00677,  34.01474,    39.948,    4.0026,  17.00737], dtype=np.float64)

# thermal conductivity [W / m / K]
lambdas = np.array([       0.0,    0.1805,       0.0,   0.02658,       0.0,     0.016,   0.02583,       0.0,    0.5863,    0.0177,     0.151,       0.0], dtype=np.float64)

index = dict(
         H= 0,     H2= 1,      O= 2,     O2= 3,     OH= 4,    H2O= 5,     N2= 6,    HO2= 7,   H2O2= 8,     AR= 9,
        HE=10,   OHEX=11
)

indexOfWater = 5
K = 12   # Number of species


"""________________________________NASA polynomials________________________________"""

N = 5    # degree of polynomials
TempRange = np.array([
    #   T_low   T_high    T_mid 
    [   200.0,  6000.0,  1000.0],    # H
    [   200.0,  6000.0,  1000.0],    # H2
    [   200.0,  6000.0,  1000.0],    # O
    [   200.0,  6000.0,  1000.0],    # O2
    [   200.0,  6000.0,  1000.0],    # OH
    [   200.0,  6000.0,  1000.0],    # H2O
    [   200.0,  6000.0,  1000.0],    # N2
    [   200.0,  5000.0,  1000.0],    # HO2
    [   200.0,  6000.0,  1000.0],    # H2O2
    [   200.0,  6000.0,  1000.0],    # AR
    [   200.0,  6000.0,  1000.0],    # HE
    [   300.0,  5000.0,  1000.0]     # OHEX
], dtype=np.float64)

# LOW NASA coefficients
a_low = np.array([
    #             a_1              a_2              a_3              a_4              a_5              a_6              a_7 
    [             2.5, -1.21434192e-16,  3.14642268e-19,   -3.308722e-22,  1.21191563e-25,        25473.66,     -0.44668285],    # H
    [      2.68763434,   0.00508352924, -1.09134795e-05,  9.75972638e-09, -2.98961778e-12,     -948.720664,    -0.706735279],    # H2
    [      3.11201699,  -0.00280464778,  5.23974549e-06, -4.42444222e-09,  1.39393321e-12,      29127.3034,      2.27964278],    # O
    [      3.66627111,  -0.00201629566,  6.94873168e-06, -6.16244023e-09,   1.7591906e-12,     -1053.52472,      4.12801194],    # O2
    [      3.98520224,  -0.00234383614,  4.44744427e-06, -3.67375964e-09,  1.27653925e-12,      3369.50653,   -0.0765438871],    # OH
    [      4.23584018,  -0.00235035878,  7.44852451e-06, -6.61473728e-09,  2.24734749e-12,     -30297.0623,    -0.999620648],    # H2O
    [      3.58256851, -0.000558781407,  7.83391479e-07,  8.73631719e-10, -7.49971666e-13,     -1051.60018,      2.75873431],    # N2
    [      3.91532945,  -0.00148787475,  1.15167354e-05, -1.25715747e-08,    4.354213e-12,      298.674804,      5.28114803],    # HO2
    [      3.87148687,   0.00289650863,  6.57196077e-06, -9.23920581e-09,  3.42064984e-12,     -17666.9584,      5.06975745],    # H2O2
    [             2.5, -1.21434192e-16,  3.14642268e-19,   -3.308722e-22,  1.21191563e-25,        -745.375,      4.37967491],    # AR
    [             2.5, -1.21434192e-16,  3.14642268e-19,   -3.308722e-22,  1.21191563e-25,        -745.375,     0.928723974],    # HE
    [       3.6372573,  0.000185157457, -1.67634008e-06,  2.38739293e-09, -8.43216802e-13,      50021.3008,      1.35889652]     # OHEX
], dtype=np.float64)

# LOW NASA coefficients
a_high = np.array([
    #             a_1              a_2              a_3              a_4              a_5              a_6              a_7 
    [             2.5,             0.0,             0.0,             0.0,             0.0,        25473.66,     -0.44668285],    # H
    [      2.93286575,  0.000826608026, -1.46402364e-07,  1.54100414e-11,   -6.888048e-16,       -816.2239,     -1.02647801],    # H2
    [      2.54363697, -2.73162486e-05,  -4.1902952e-09,  4.95481845e-12, -4.79553694e-16,      29226.5295,      4.92264671],    # O
    [      3.66096065,  0.000656365811, -1.41149627e-07,  2.05797935e-11, -1.29913436e-15,     -1214.90829,      3.41609021],    # O2
    [      2.83853033,   0.00110741289, -2.94000209e-07,  4.20698729e-11,  -2.4228989e-15,      3697.87047,      5.84498899],    # OH
    [       2.6770389,    0.0029731816,  -7.7376889e-07,   9.4433514e-11,  -4.2689991e-15,     -29886.2362,      6.88231731],    # H2O
    [      2.95257637,    0.0013969004, -4.92631603e-07,  7.86010195e-11, -4.60755204e-15,     -924.423063,      5.87156477],    # N2
    [      4.17228741,   0.00188117627, -3.46277286e-07,  1.94657549e-11,  1.76256905e-16,      34.5761323,      2.96009637],    # HO2
    [      4.57977305,   0.00405326003,  -1.2984473e-06,    1.982114e-10, -1.13968792e-14,     -18003.0958,      0.66774843],    # H2O2
    [             2.5,             0.0,             0.0,             0.0,             0.0,        -745.375,      4.37967491],    # AR
    [             2.5,             0.0,             0.0,             0.0,             0.0,        -745.375,     0.928723974],    # HE
    [         2.88273,    0.0010139743,   -2.276877e-07,    2.174683e-11,   -5.126305e-16,      50301.4063,      5.59571607]     # OHEX
], dtype=np.float64)


"""________________________________Reaction constants________________________________"""

I = 30    # Number of reactions
# Pre-exponential factors [cm^3/mol/s v 1/s]
A = np.array([
      5071200000000000.0,           1255400.0,          13193000.0,             84999.0,          4.9806e+18,
      6165000000000000.0,           4.714e+18,          1.4818e+24,     4650000000000.0,           2123100.0,
        57734000000000.0,    32500000000000.0,      958400000000.0,      130000000000.0,  1604800000000000.0,
                214800.0,    24100000000000.0,          9.7543e+19,           9550000.0,     1740000000000.0,
        75900000000000.0,    15000000000000.0,     5930000000000.0,     2950000000000.0,      108000000000.0,
         6010000000000.0,     1310000000000.0,     1690000000000.0,           1450000.0,     2100000000000.0
], dtype=np.float64)

# Temperature exponent [-]
b = np.array([
                -0.48596,             2.27039,             1.87803,             2.26419,            -1.21273,
                    -0.5,                -1.0,            -2.53792,                0.44,              2.1133,
                     0.0,                 0.0,             0.42008,                 0.0,                 0.0,
                  2.3219,                 0.0,            -1.92495,                 2.0,                 0.0,
                     0.0,                 0.0,                 0.5,                 0.5,                 0.5,
                     0.5,                 0.5,                 0.0,                 0.0,                 0.5
], dtype=np.float64)

# Activation energy [cal/mol]
E = np.array([
             16128.31392,          6957.58464,          3151.30176,         -1784.96266,           612.09734,
                     0.0,                 0.0,           120.79791,                 0.0,          -1624.8937,
                171.0383,                 0.0,          -948.68928,         -1630.15978,         15551.03232,
             -3402.70243,          3970.40573,          9426.08448,          3970.40573,           318.03149,
              7269.73402,          5975.60976,          -860.08003,          -444.03984,         -1242.11923,
               -764.0784,          -167.02416,          4135.42282,                 0.0,          -478.04083
], dtype=np.float64)


"""________________________________Reaction matrixes________________________________"""

# Forward reaction matrix
nu_forward = np.array([
    #   H   H2    O   O2   OH  H2O   N2  HO2 H2O2   AR   HE OHEX 
    [   1,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0],    #  0. H+O2=O+OH
    [   0,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0],    #  1. O+H2=H+OH
    [   0,   1,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0],    #  2. OH+H2=H+H2O
    [   0,   0,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0],    #  3. 2OH=O+H2O
    [   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0],    #  4. 2H+M=H2+M
    [   0,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0],    #  5. 2O+M=O2+M
    [   1,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0],    #  6. O+H+M=OH+M
    [   1,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0],    #  7. H+OH+M=H2O+M
    [   1,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0],    #  8. H+O2(+M)=HO2(+M)
    [   1,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0],    #  9. H+HO2=H2+O2
    [   1,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0],    # 10. HO2+H=2OH
    [   0,   0,   1,   0,   0,   0,   0,   1,   0,   0,   0,   0],    # 11. HO2+O=OH+O2
    [   0,   0,   0,   0,   1,   0,   0,   1,   0,   0,   0,   0],    # 12. HO2+OH=H2O+O2
    [   0,   0,   0,   0,   0,   0,   0,   2,   0,   0,   0,   0],    # 13. 2HO2=H2O2+O2
    [   0,   0,   0,   0,   0,   0,   0,   2,   0,   0,   0,   0],    # 14. 2HO2=H2O2+O2
    [   0,   0,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0],    # 15. 2OH(+M)=H2O2(+M)
    [   1,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0],    # 16. H2O2+H=H2O+OH
    [   1,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0],    # 17. H2O2+H=H2+HO2
    [   0,   0,   1,   0,   0,   0,   0,   0,   1,   0,   0,   0],    # 18. H2O2+O=OH+HO2
    [   0,   0,   0,   0,   1,   0,   0,   0,   1,   0,   0,   0],    # 19. H2O2+OH=H2O+HO2
    [   0,   0,   0,   0,   1,   0,   0,   0,   1,   0,   0,   0],    # 20. H2O2+OH=H2O+HO2
    [   1,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0],    # 21. H+O+M=OHEX+M
    [   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   1],    # 22. OHEX+H2O=OH+H2O
    [   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1],    # 23. OHEX+H2=OH+H2
    [   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   1],    # 24. OHEX+N2=OH+N2
    [   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   1],    # 25. OHEX+OH=2OH
    [   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1],    # 26. OHEX+H=OH+H
    [   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   1],    # 27. OHEX+AR=OH+AR
    [   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1],    # 28. OHEX=OH+HV
    [   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   1]     # 29. OHEX+O2=OH+O2
], dtype=np.float64)

# Backward reaction matrix
nu_backward = np.array([
    #   H   H2    O   O2   OH  H2O   N2  HO2 H2O2   AR   HE OHEX 
    [   0,   0,   1,   0,   1,   0,   0,   0,   0,   0,   0,   0],    #  0. H+O2=O+OH
    [   1,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0],    #  1. O+H2=H+OH
    [   1,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0],    #  2. OH+H2=H+H2O
    [   0,   0,   1,   0,   0,   1,   0,   0,   0,   0,   0,   0],    #  3. 2OH=O+H2O
    [   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0],    #  4. 2H+M=H2+M
    [   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0],    #  5. 2O+M=O2+M
    [   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0],    #  6. O+H+M=OH+M
    [   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0],    #  7. H+OH+M=H2O+M
    [   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0],    #  8. H+O2(+M)=HO2(+M)
    [   0,   1,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0],    #  9. H+HO2=H2+O2
    [   0,   0,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0],    # 10. HO2+H=2OH
    [   0,   0,   0,   1,   1,   0,   0,   0,   0,   0,   0,   0],    # 11. HO2+O=OH+O2
    [   0,   0,   0,   1,   0,   1,   0,   0,   0,   0,   0,   0],    # 12. HO2+OH=H2O+O2
    [   0,   0,   0,   1,   0,   0,   0,   0,   1,   0,   0,   0],    # 13. 2HO2=H2O2+O2
    [   0,   0,   0,   1,   0,   0,   0,   0,   1,   0,   0,   0],    # 14. 2HO2=H2O2+O2
    [   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0],    # 15. 2OH(+M)=H2O2(+M)
    [   0,   0,   0,   0,   1,   1,   0,   0,   0,   0,   0,   0],    # 16. H2O2+H=H2O+OH
    [   0,   1,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0],    # 17. H2O2+H=H2+HO2
    [   0,   0,   0,   0,   1,   0,   0,   1,   0,   0,   0,   0],    # 18. H2O2+O=OH+HO2
    [   0,   0,   0,   0,   0,   1,   0,   1,   0,   0,   0,   0],    # 19. H2O2+OH=H2O+HO2
    [   0,   0,   0,   0,   0,   1,   0,   1,   0,   0,   0,   0],    # 20. H2O2+OH=H2O+HO2
    [   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1],    # 21. H+O+M=OHEX+M
    [   0,   0,   0,   0,   1,   1,   0,   0,   0,   0,   0,   0],    # 22. OHEX+H2O=OH+H2O
    [   0,   1,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0],    # 23. OHEX+H2=OH+H2
    [   0,   0,   0,   0,   1,   0,   1,   0,   0,   0,   0,   0],    # 24. OHEX+N2=OH+N2
    [   0,   0,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0],    # 25. OHEX+OH=2OH
    [   1,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0],    # 26. OHEX+H=OH+H
    [   0,   0,   0,   0,   1,   0,   0,   0,   0,   1,   0,   0],    # 27. OHEX+AR=OH+AR
    [   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0],    # 28. OHEX=OH+HV
    [   0,   0,   0,   1,   1,   0,   0,   0,   0,   0,   0,   0]     # 29. OHEX+O2=OH+O2
], dtype=np.float64)

nu = nu_backward - nu_forward

reaction_order = np.array([
    [   2],    #  0. H+O2=O+OH
    [   2],    #  1. O+H2=H+OH
    [   2],    #  2. OH+H2=H+H2O
    [   2],    #  3. 2OH=O+H2O
    [   2],    #  4. 2H+M=H2+M
    [   2],    #  5. 2O+M=O2+M
    [   2],    #  6. O+H+M=OH+M
    [   2],    #  7. H+OH+M=H2O+M
    [   2],    #  8. H+O2(+M)=HO2(+M)
    [   2],    #  9. H+HO2=H2+O2
    [   2],    # 10. HO2+H=2OH
    [   2],    # 11. HO2+O=OH+O2
    [   2],    # 12. HO2+OH=H2O+O2
    [   2],    # 13. 2HO2=H2O2+O2
    [   2],    # 14. 2HO2=H2O2+O2
    [   2],    # 15. 2OH(+M)=H2O2(+M)
    [   2],    # 16. H2O2+H=H2O+OH
    [   2],    # 17. H2O2+H=H2+HO2
    [   2],    # 18. H2O2+O=OH+HO2
    [   2],    # 19. H2O2+OH=H2O+HO2
    [   2],    # 20. H2O2+OH=H2O+HO2
    [   2],    # 21. H+O+M=OHEX+M
    [   2],    # 22. OHEX+H2O=OH+H2O
    [   2],    # 23. OHEX+H2=OH+H2
    [   2],    # 24. OHEX+N2=OH+N2
    [   2],    # 25. OHEX+OH=2OH
    [   2],    # 26. OHEX+H=OH+H
    [   2],    # 27. OHEX+AR=OH+AR
    [   1],    # 28. OHEX=OH+HV
    [   2]     # 29. OHEX+O2=OH+O2
], dtype=np.int64)


"""________________________________Three-body reactions________________________________"""

ThirdBodyIndexes = np.array([   4,   5,   6,   7,   8,  15,  21], dtype=np.int64)
ThirdBodyCount = 7

# third-body efficiency factors
alfa = np.array([
    #       H       H2        O       O2       OH      H2O       N2      HO2     H2O2       AR       HE     OHEX 
    [     1.0,     2.5,     1.0,     1.0,     1.0,    12.0,     1.0,     1.0,     1.0,     1.0,    0.83,     1.0],    #  4. 2H+M=H2+M
    [     1.0,     2.5,     1.0,     1.0,     1.0,    12.0,     1.0,     1.0,     1.0,    0.83,    0.83,     1.0],    #  5. 2O+M=O2+M
    [     1.0,     2.5,     1.0,     1.0,     1.0,    12.0,     1.0,     1.0,     1.0,    0.75,    0.75,     1.0],    #  6. O+H+M=OH+M
    [     1.0,     2.5,     1.0,     1.0,     1.0,    12.0,     1.0,     1.0,     1.0,    0.38,    0.44,     1.0],    #  7. H+OH+M=H2O+M
    [     1.0,   1.511,     1.0,     1.0,     1.0,  11.372,     1.0,     1.0,     1.0,   0.474,    0.65,     1.0],    #  8. H+O2(+M)=HO2(+M)
    [     1.0,    2.47,     1.0,     0.8,     1.0,     5.0,     1.0,     1.0,    5.13,    0.67,    0.43,     1.0],    # 15. 2OH(+M)=H2O2(+M)
    [     1.0,     1.0,     1.0,     0.4,     1.0,     6.5,     0.4,     1.0,     1.0,    0.35,     1.0,     1.0]     # 21. H+O+M=OHEX+M
], dtype=np.float64)


"""________________________________Irreversible reactions________________________________"""

IrreversibleIndexes = np.array([], dtype=np.int64)
IrreversibleCount = 0


"""________________________________Pressure-dependent reactions________________________________"""

PressureDependentIndexes = np.array([   8,  15], dtype=np.int64)
PressureDependentCount = 2

LindemannIndexes = np.array([], dtype=np.int64)
LindemannCount = 0

# Fall-off parameters
ReacConst = np.array([
    #               A_0                b_0                E_0 
    [        5.2669e+19,          -1.37367,               0.0],    #  8. H+O2(+M)=HO2(+M)
    [        1.9928e+18,          -1.17797,      -4273.096032]     # 15. 2OH(+M)=H2O2(+M)
], dtype=np.float64)

TroeIndexes = np.array([   8,  15], dtype=np.int64)
TroeCount = 2

# Troe parameters
Troe = np.array([
    #              alfa               T***                 T*                T** 
    [              0.67,             1e-30,             1e+30,             1e+30],    #  8. H+O2(+M)=HO2(+M)
    [              0.43,             1e-30,             1e+30,             1e+30]     # 15. 2OH(+M)=H2O2(+M)
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

