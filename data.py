"""
This python file contains numerical data missing from the .inp OpenSmoke files.
Recommended usage:
   importing: import data
   usage: print('Molar mass of hydrogen: =', data.W['H'])
"""

# thermal conductivity [W/m/K]
lambdas = dict(
    HE=0.1513,
    NE=0.0491,
    AR=0.01772,
    KR=0.00943,
    H2=0.1805,
    O2=0.02658,
    H2O=0.598,
    H2O2=0.5863,
    O3=0.019854,
    CO=0.024,
    CO2=0.01663,
    N2=0.02583,
    NH3=0.00244,
    NO2=0.00988,
)

header='''\"\"\"
This is an automatically generated file.
This python file contains all the numerical constants and coefficents required for the simulations.
Recommended usage:
   importing: import parameters as par
   usage: print('Universal gas constant =', par.R_g)
\"\"\"'''

physical_constants = '''c_L = 1483.0       # Liquid sound speed [m/s]
rho_L = 998.2      # Liquid density [kg/m^3]
sigma = 71.97e-3   # Surface tension [N/m]
mu_L = 0.001       # Dynamic viscosity [Pa*s]
R_v = 461.5227     # Specific gas constant of water [J/kg/K]

R_g = 8.31446      # Universal gas constant [J/mol/K]
R_erg = 8.31446e7  # Universal gas constant [erg/mol/K]
R_cal = 1.987      # Universal gas constant [cal/mol/K]
N_A = 6.02214e23   # Avogadro's number [-]
h = 6.62607015e-34 # Planck constant [m^2*kg/s]'''

valid_elements='''H, HE, LI, BE, B, C, N, O, F, NE, NA, MG, AL, SI, P, S, CL, AR, K, CA, SC, TI, V, CR, MN, FE, CO, NI, CU, ZN, GA,
GE, AS, SE, BR, KR, RB, SR, Y, ZR, NB, MO, TC, RU, RH, PD, AG, CD, IN, SN, SB, TE, I, XE, CS, BA, LA, CE, PR, ND,
PM, SM, EU, GD, TB, DY, HO, ER, TM, YB, LU, HF, TA, W, RE, OS, IR, PT, AU, HG, TL, PB, BI, PO, AT, RN, FR, RA,
AC, TH, PA, U, NP, PU, AM, CM, BK, CF, ES, FM, D, E'''

# molar masses [g/mol]
W = dict(
  # 1. row
    H=1.00797,
    HE=4.00260,
  # 2. row   
    LI=6.941,
    BE=9.01218,
    B=10.81,
    C=12.011,
    N=14.0067,
    O=15.9994,
    F=18.998403,
    NE=20.179,
  # 3.row
    NA=22.98977,
    MG=24.305,
    AL=26.98154,
    SI=28.0855,
    P=30.97376,
    S=32.06,
    CL=35.453,
    AR=39.948,
  # 4. row
    K=39.0983,
    CA=40.08,
    SC=44.9559,
    TI=47.90,
    V=50.9415,
    CR=51.996,
    MN=54.9380,
    FE=55.847,
    CO=58.9332,
    NI=58.70,
    CU=63.546,
    ZN=65.38,
    GA=69.72,
    GE=72.59,
    AS=74.9216,
    BR=79.904,
    KR=83.80,
  # 5. row
    RB=85.4678,
    SR=87.62,
    Y=88.9059,
    ZR=91.22,
    NB=92.9064,
    MO=95.94,
    TC=98.0,
    RU=101.07,
    RH=102.9055,
    PD=106.4,
    AG=107.868,
    CD=112.41,
    IN=114.82,
    SN=118.69,
    SB=121.75,
    TE=127.60,
    I=126.9045,
    XE=131.30,
  # 6. row
    CS=132.9054,
    BA=137.33,
    LA=138.9055,
    CE=140.12,
    PR=140.9077,
    ND=144.24,
    PM=145.0,
    SM=150.4,
    EU=151.96,
    GD=157.25,
    TB=158.9254,
    DY=162.50,
    HO=164.9304,
    ER=167.26,
    TM=168.9342,
    YB=173.04,
    LU=174.967,
    HF=178.49,
    TA=180.9479,
    W=183.85,
    RE=186.207,
    OS=190.2,
    IR=192.22,
    PT=195.09,
    AU=196.9665,
    HG=200.59,
    TL=204.37,
    PB=207.2,
    BI=208.9804,
    PO=209.0,
    AT=210.0,
    RN=222.0,
  # 7. row
    FR=223.0,
    RA=226.0254,
    AC=227.0278,
    TH=232.0381,
    PA=231.0359,
    U=238.029,
    NP=237.0482,
    PU=242.0,
    AM=243.0,
    CM=247.0,
    BK=247.0,
    CF=251.0,
    ES=252.0,
    FM=257.0,
  # Other
    D=2.014, # deuterium
    E=5.4858e-4, # electron
)