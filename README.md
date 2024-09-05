#  Bubble dynamics simulation

Tutorial video: [Magyar](https://youtu.be/YWsT1ktUzVw), [English]()

About the HDS Sonochemistry research group:
 * University: Budapest University of Technology and Economics
 * Faculty of Mechanical Engineering
 * Department of Hydrodynamic Systems
 * [web page](https://www.hds.bme.hu/research/BubbleDynamics/index.html) (hun)

About the author:
 * name: Áron Kozák
 * email: kozi0223@gmail.com
 * degrees:
 	* Mechatronics Engineering BSc
	* Mechanical Engineering Modelling MSc in progress

## Table of content

1. [**Python basics**](#python_basics): Feel free to skip this part. Helps with installations, briefly introduces some special python construct, and links tutorials. Might be helpful, if you are new to python.

2. [**The simulation**](#the_simulation): Introduces some of the basic functions in [full_bubble_model.py](full_bubble_model.py). Shows, how to solve the differential equation with different control parameters. Shows how to plot, process, or save the results.

3. [**Bruteforce parameter sweep**](#bruteforce_parameter_sweep): A simple multithread example. A control parameter space is sweeped with the given resolution. Results are saved into CSV files.

4. [**Read CSV files**](#read_csv_files): This example showcases how to read CSV files from a given folder. The data is loaded into a pandas dataframe. Different statistics and data manipulation can be done with the dataframe.

5. [**Gradient descent**](#gradient_descent): This 2 examples showcases a global optimazition algorithm. It is way faster than a bruteforce parameter sweep, exspecially with a higher number of dimensions. The search can be visualized on a 2D contour plot.

6. [**Create plots**](#create_plots): This example is about creating simple plots using the simulation.

<a name="python_basics"/>

##  Python basics

Feel free to skip this part. If you are absolutely new to python, you may start by checking out [Learn Python](https://www.learnpython.org/) or [w3schools](https://www.w3schools.com/python/) or watch a [1 hour video](https://www.youtube.com/watch?v=kqtD5dpn9C8) or a slightly [longer video](https://www.youtube.com/watch?v=eWRfhZUzrAc) with more details

### Installing python

First download [Python](https://www.python.org/downloads/), and make sure you, that you install a version compatible with [numba](https://numba.readthedocs.io/en/stable/user/installing.html). Don't forget to check the *Add python.exe to PATH* box in the installation dialogue. Open a command prompt, and try the `python --version` command.

Then check that [pip](https://pip.pypa.io/en/stable/installation/) was automatically installed with `pip --version`. If not, you may use `python -m ensurepip --upgrade` command. Pip is used to install packages/libraries. You can run pip commands in a python notebook's code cell with an exclamation mark (`!pip...`). Basic commands:
* install: `pip install <package_name>`
* uninstall: `pip uninstall <package_name>`
* upgrade version: `pip install <package_name> --upgrade`
* check version: `pip show <package_name>`

You may start by the following installations:
~~~
pip install numpy matplotlib scipy numba func_timeout pandas termcolor
~~~

Finally you need to choose one of the several IDEs available. The most popular choice is [Visual Studio Code](https://code.visualstudio.com/download). An alternative is Jupyter Lab. You can downlad the new [desktop version](https://jupyter.org/), or download [Anaconda](https://www.anaconda.com/products/distribution), and launch it from the Anaconda Prompt by typing `jupyter lab`. This way only one drive will be avalible. Launch with `jupyter lab --notebook-dir=F:/` to change to *F* drive.

Note that there are 2 types of python file out there:
* **.py** files: The original python format. It's just plain text, like *.c* *.cpp* *.h*. You can edit them in any text editor, and run them from a command prompt without any IDEs. Only *.py* files can be imported into other projects as libraries. **The .py files contain the source code of the package.**
* **.ipynb** files: The notebook format, in a text editor you will see a html like structure, which has to be rendered, thus you need to use an IDE. A notebook contains code cells, which can run separately, and have they own output. There are also markdown cells, which displays formatted text with headers, latex formulas, images, etcetera. Some IDEs can generate a table of content from the headers. **The .ipynb files contain examples to showcase the capabilities of the package.**

### The f-string

You can use f-strings to easily format strings:
~~~Python
a=4; pi=3.14159265128

# print variables:
print(f'a={a} or {a=}')

# format the output:
print(f'unformatted: {pi}')
print(f'set decimal places: {pi: .4f}')
print(f'set width: |{pi: <8.4}| or |{pi: >8.4}|')

# format using a variable as number of decimal places:
print(f'set decimal places: {pi: .{a}f}')
~~~

~~~
a=4 or a=4
unformatted: 3.14159265128
set decimal places:  3.1416
set width: |3.142   | or |   3.142|
set decimal places:  3.1416
~~~

### About for loops

Python have some very convenient syntax for loops:
~~~Python
a = ['zero', 'one', 'two']
b = ['nulla', 'egy', 'ketto']

# regular for loop:
for i in range(len(a)):
    print(a[i], end=' ')
    
# iterate trought elements
print(' ')
for element in a:
    print(element, end=' ') # by modifying element, a wont change
    
# use enumerate
print('\n')
for i, element in enumerate(a):
    print(i, element) # modifying a[i] will change a, but modifying element won't
    
# use zip
print('')
for element1, element2 in zip(a, b):
    print(element1, element2)
~~~

~~~
zero one two  
zero one two 

0 zero
1 one
2 two

zero nulla
one egy
two ketto
~~~

List comprehensions:
~~~Python
# create a list with a for loop:
a = []
for digit in range(10):
    a.append(digit)
print(a)

# or use a shortcut:
a = [digit for digit in range(10)]
print(a)

# get the square of all even digits:
a = [digit**2 for digit in range(10) if digit%2==0]
print(a)
~~~

~~~
[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
[0, 4, 16, 36, 64]
~~~

###  The dictionary

Dictionaries are data containers used to store different types of variables in an organized manner. Each member variable can be accessed via a key. It's similar to STL maps in c++. Basic usage:

~~~Python
# create a dictionary
datas = dict(
	data1=1,
	data2='one',
)

# add new element
datas['data3']=1.0

# remove element
datas.pop('data1')

# print dictionary
print(datas)

# iterate trought key
for key in datas:
	print(datas[key], end=' ')
	
# acces element
print(datas['data2'])

# check if a key is in the dictionary
print('data2' in datas, 'data1' in datas)
~~~

~~~
{'data2': 'one', 'data3': 1.0}
one 1.0
'one'
True False
~~~

For easier syntax, you can use the dot notation to acces an element via this 4 lines of magic:

~~~Python
class dotdict(dict):
	__getattr__ = dict.get
	__setattr__ = dict.__setitem__
	__delattr__ = dict.__delitem__

datas = dotdict(dict( data=1 ))
datas.data
~~~

~~~
1
~~~

This code is also included in *full_bubble_mode.py*. Other features of the dictionary is not affected. Most of my codes uses the dot notation. Note, that this makes the dictionaries a bit slower. Also, multiprocessing.Pool() doesn't like dotdicts.

<a name="the_simulation"/>

##  The simulation

### Parameters, automatic generation

Before running any simulations, the numerical constants must be extracted. They are provided by the ELTE Chemical Kinetics Laboratory in the form of OpenSmoke files with .inp file extension. These are text files with a barely human readable and slightly inconsistent formatting. Lets use [chem_Otomo2018_without_O.inp](<INP file examples/chem_Otomo2018_without_O.inp>), which contains 12 nitrogen and hydrogen compounds with 35 reactions. You may trun this, or any other reaction model into a python file using the `extract()` function from [inp_data_extractor.py](inp_data_extractor.py):

~~~Python
from Bubble_dynamics_simulation import inp_data_extractor as inp
path = 'Bubble_dynamics_simulation/INP file examples/chem_Otomo2018_without_O.inp'
inp.extract(path)
~~~

~~~
path=Bubble_dynamics_simulation/INP file examples/chem_Otomo2018_without_O.inp
Note, lambda value for specie 'H' is not in data.py: 0.0 is used
...
Warning, third body 'AR' is not in species in line 64 (' H2/2.5/ H2O/12/ AR/0.0/ ') in reaction 'H2+M=H+H+M'
model: chem_Otomo2018_without_O
File 'parameters.py' succesfully created
~~~

As a result, a new file, [parameters.py](parameters.py), is generated in the current working directory. Error messages and warnings may also be displayed. Some constants are not includ-ed in .inp files, they come from [data.py](data.py). The generated file contains all the required constants in an organized manner. To check it, open in a text editor or run:

~~~Python
import parameters as par
print(par.model)
print(par.species)
~~~

~~~
chem_Otomo2018_without_O
['NH3' 'H2' 'H' 'NH2' 'NH' 'N' 'NNH' 'N2H4' 'N2H3' 'N2H2' 'H2NN' 'N2']
~~~

### Control parameters

Now, the simulation, [full_bubble_model.py](full_bubble_model.py) can be loaded. Beforehand, set the varia-bles in this file under the `___settings___` comment. These settings are displayed once the file is imported, and cannot be changed later. For this step, having parameters.py is required.

~~~Python
from Bubble_dynamics_simulation import full_bubble_model as de
~~~

~~~
model: chem_Otomo2018_without_O
target specie: NH3
excitation: sin_impulse (control parameters: ['p_A', 'freq', 'n'])
enable heat transfer: True enable evaporation: False enable reactions: True enable dissipated energy: True
~~~

The control parameters, which can be used to influence the system, are passed to the simulation in a dictionary. These control parameters can be tuned to lower the energy demand of ammonia production. To obtain such a dictionary you may use `cpar = de.example_cpar()`. You may print this dictionary and use it as code with `de.print_cpar(de.example_cpar())`. The result looks something like this:

~~~Python
T = 293.15 # [K]
cpar = dict(
    ID = 0,                                         # ID of control parameter
  # Initial conditions:
    R_E = 100.0e-6,                                 # bubble equilibrium radius [m]
    ratio = 1.00,                                   # R_0/R_E [-]
    gases = [par.index['H2'], par.index['N2']],     # indexes of species
    fractions = [0.75, 0.25],                       # molar fractions of species
  # Ambient parameters:
    P_amb = 101325.00,                              # ambient pressure [Pa]
    T_inf = T,                                      # ambient temperature [K]
  # Liquid parameters:
    alfa_M = 0.35,                                  # water accommodation coefficient [-]
    P_v = de.vapour_pressure(T),                     # vapour pressure [Pa]
    mu_L = de.viscosity(T),                         # dynamic viscosity [Pa*s]
    rho_L =  998.20,                                # density [kg/m^3]
    c_L = 1483.00,                                  # sound speed [m/s]
    surfactant = 1.00,                              # surface tension modfier [-]
  # Excitation parameters:
    p_A = -2.0e5,                                   # [Pa]
    freq = 10000.0,                                 # [Hz]
    n = 1.0,                                        # [-]
)
~~~

Note, that some control parameters, `P_v` and `mu_L` are calculated from the temperature. To access a data member, use indexing: `cpar['R_E']`. 

### Numerical solution

In order to obtain the numerical solution, use the `solve()` function. The only mandatory argument is `cpar`; however, additional settings are available, such as the simulated time interval, and timeout limits for the LSODA and Radau solvers. The function checks the validity of cpar, calculates the initial conditions, and run the mentioned numeric solvers with error management and runtime measurements:

~~~Python
num_sol, error_code, elapsed_time = de.solve(cpar)
~~~

To access the numerical solution, use `num_sol.t` and `num_sol.y`. To save the results, check out the `Make_dir` class. A post processing function called `get_data()` is also available. The function returns another dictionary containing the control parameters, errors, computational time, length of the numeric solution, the initial conditions, the final time step, the energy dissipated or put into the system, and the energy demand of ammonia along with several additional information, like the peak temperature. Usually, only this dictionary is saved in a CSV (comma separated values) text file. You may call post processing as:

~~~Python
data = de.get_data(cpar, num_sol, error_code, elapsed_time)
~~~

To format and display this data, use `de.print_data(cpar, data)`. 

### Plot simulations

There is a built-in function, that calculates the initial conditions, solves the ODE, plots the results, and displays the post processing data. The name of the function is `plot()`, and while it has several optional arguments, the only mandatory input is cpar:

~~~Python
de.plot(cpar)
~~~

~~~
succecfully solved with LSODA solver
~~~
![Image not found](https://github.com/hihihi2001/Bubble_dynamics_simulation/assets/42745647/a80ba305-85cf-4500-a936-474bfcd05b40)
![Image not found](https://github.com/hihihi2001/Bubble_dynamics_simulation/assets/42745647/9f6281e0-17e6-4941-810b-bbb4fb4e753c)
~~~
Control parameters:
    ID = 0,                                      # ID of control parameter (not used during calculation)
  # Initial conditions:
    R_E =  0.00010000,                           # bubble equilibrium radius [m]
    ratio =  1.00,                               # initial radius / equilibrium radius R_0/R_E [-]
    gases = [par.index['H2'], par.index['N2']],  # indexes of species in initial bubble (list of species indexes)
    fractions = [0.75, 0.25],                    # molar fractions of species in initial bubble (list of fractions for every gas)
  # Ambient parameters:
    P_amb =  101325.00,                          # ambient pressure [Pa]
    T_inf =  293.15,                             # ambient temperature [K]
  # Liquid parameters:
    alfa_M =  0.3500,                            # water accommodation coefficient [-]
    P_v =  2338.34,                              # vapour pressure [Pa]
    mu_L =  0.0010,                              # dynamic viscosity [Pa*s]
    rho_L =  998.20,                             # density [kg/m^3]
    c_L =  1483.00,                              # sound speed [m/s]
    surfactant =  1.00,                          # surfactant (surface tension modfier) [-]
  # Excitation parameters: (excitation_type = sin_impulse)
    p_A = -200000.00,                            # [Pa]
    freq =  10000.00,                            # [Hz]
    n =  1.00,                                   # [-]

Simulation info:
    error_code = 0 (success = True)
    elapsed_time = 0.51 [s]
    steps = 13854 [-]
Final state:
    R_final = 99.70 [um];   R_dot_final =-4.212719033050712e-12 [m/s];   T_final = 293.15 [K]
    n_NH3_final = 1.56e-12 [mol]
    Final molar concentrations: [mol/cm^3]
        NH3   :  3.753749e-07;    H2    :  3.134108e-05;    H     : -4.106202e-14;    NH2   :  1.564072e-15;    
        NH    : -1.583905e-24;    N     :  5.789647e-12;    NNH   :  7.193518e-22;    N2H4  :  3.105080e-12;    
        N2H3  :  6.187001e-11;    N2H2  : -7.644744e-13;    H2NN  : -2.723146e-14;    N2    :  1.044699e-05;    
        
Results:
    collapse_time = 8.28390263916771e-05 [s]
    T_max = 5129.74 [K]
    expansion work = 0.0 [J]
    dissipated acoustic energy = 3.221739162989513e-05 [J]
    energy demand = 1213.8764401269318 [MJ/kg of NH3]
~~~

### Save simulation results

Use the `Make_dir` class to save results into CSV files. To use, first run:

~~~ Python
file = de.Make_dir('test')
~~~

This will create a folder called test, or use the existing one. Inside the folder, you can create a new, automatically numbered CSV file with a header:

~~~ Python
file.new_file()
~~~

The default name is *output_1.csv* with the number automatically incrementing itself. You can add a new line and store results:

~~~Python
data = de.simulate(cpar)
file.write_line(data)
~~~

The `simulate()` functions automatically runs `solve()` and `get_data()`, and only returns `data`. It uses regular dictionaries, thus it can be used with the multiprocessing module. You can close the CSV to save it:

~~~Python
file.close()
~~~

You can also save the numerical solution:

~~~Python
num_sol, error_code, elapsed_time = de.solve(cpar)
data = de.get_data(cpar, num_sol, error_code, elapsed_time)
file.write_solution(data, num_sol, 'testfile')
~~~

Two files will be created in the *test* folder. *testfile_data.csv* contains the data dictionary, while *testfile_num_sol.csv* contains the numerical solution. This function is independent of `new_file()` and `close()`.

You may read the detailed description of any public function (not starting with an underscore) in the source code or in [DOCUMENTATION.md](DOCUMENTATION.md).

<a name="bruteforce_parameter_sweep"/>

##  Bruteforce parameter sweep

Example file: [Bruteforce_parameter_sweep.ipynb](Bruteforce_parameter_sweep.ipynb)

|![Image not found](https://github.com/hihihi2001/Bubble_dynamics_simulation/assets/42745647/db01c983-1504-4e3b-9d35-f74aaa71d716)|
|:---:|
|Visualization of a bruteforce parameter sweep in 2 dimensions. The units are meter and Pascal respectively. The simulation is evaluated in all 400 balck points, and the one with the smallest energy demand is choosen. |

In this example, you can create a parameter study. You give a list for all control parameters, containing all their possible values. Then the program makes a long list of all the possible combinations. Affterward importing [multiprocessing](https://docs.python.org/3/library/multiprocessing.html#using-a-pool-of-workers), you can solve the simulation for all combinations using all your computer cores. (e.g. R_E=[1, 2, 3, 4] [um]; expansion ratio=[5, 10, 15, 20] [-] --> total 4*4=16 combinations)

For a very simple example of multiprocessing, create a new file, called *add.py*. Inside create a function, that adds the first 2 elements of a list, but very slowly:

~~~Python
def add(inp):
    time.sleep(1)   # make it slow
    return inp[0]+inp[1]
~~~

Then in a new notebook copy this:

~~~Python
from add import add
from multiprocessing import Pool

# possible parameter values:
range1 = [1, 2, 3, 4]
range2 = [1, 2, 3, 4]

# all combinations:
combinations = []
for par1 in range1:
    for par2 in range2:
        combinations.append([par1, par2])
        
# run all combinations using multiprocessing instead of a for loop:
sums = []
with Pool() as pool:
    results = pool.imap(add, combinations)
    
    for result in results:   # process the results
        sums.append(result)
        
sums
~~~

~~~
[2, 3, 4, 5, 3, 4, 5, 6, 4, 5, 6, 7, 5, 6, 7, 8]
~~~

Normally, using a for loop, it would take 16 seconds to run all the combinations, but using all your cpu threads, you can get it done much faster. The example might seem a bit more complicated, but the core mechanic is the same. The examlpe uses pool.imap_unordered, instead of pool.imap to get result faster, but in arbitrary order. The example uses the `Make_dir` class to save the results to .csv files.

<a name="read_csv_files"/>

##  Read CSV files

Example file: [Read_CSV_files.ipynb](Read_CSV_files.ipynb)

After running a parameter study, you will have several CSV files with the results. This example uses [pandas](https://pandas.pydata.org/docs/), a popular data science library. The program lists and loads all CSV files from the given folder into a dataframe, and makes some data manipulations, like filtering solutions, and also gets some statistics. You can also retrive `cpar` and plot interesting simulations.

You can create a pandas dataframe from a CSV: `df=pandas.read_csv('example.csv')`, or from a dictionary using `df=pandas.DataFrame(example_dict)`:

~~~Python
import pandas as pd

df = pd.DataFrame(dict(num=[1, 2], name=['one', 'two']))
print(df)
~~~

~~~
   num name
0    1  one
1    2  two
~~~

The most basic data manipulation is adding new rows or columns:

~~~Python
# add new row:
new_row = pd.DataFrame(dict(num=[3], name=['three']))
df = pd.concat([df, new_row]).reset_index(drop=True)

# create a new column:
df['square'] = df['num']**2

print(df)
~~~

~~~
   num   name  square
0    1    one       1
1    2    two       4
2    3  three       9
~~~

An other very important tool is filtering a dataframe. Here we don't like the number two. To filter by multiple conditions, use bitwise and (&) / or (|):

~~~Python
print(df.loc[ (df['num']!=2) & (df['name']!='two') ])
~~~

~~~
   num   name  square
0    1    one       1
2    3  three       9
~~~

<a name="gradient_descent"/>

##  Gradient descent

|![Image not found](https://github.com/hihihi2001/Bubble_dynamics_simulation/assets/42745647/3e4c7885-bc66-4663-b777-c464ae9be620)|
|:---:|
|2D visualization of the gradient descent based optimization. The units are meter and Pascal respectively. |

Example files: [Gradient_descent.ipynb](Gradient_descent.ipynb); [Gradient_descent_2D.ipynb](Gradient_descent_2D.ipynb)

### Gradient_descent.ipynb

In this example you can use a gradient descent based algorithm to more efficiently optimize the control parameters to minimize energy demand. Basically, the program can start in any point and make steps always opposite to the local gradient. The gradient is approxiamted with finite differences. The gradient is normed, and a step size is controled seperately.

The algorithm can minimize any number that is provided in tha data dictionary returned by the `get_data()` function. You have to set the key via the `to_optimize` variable:

~~~Python
save_path = 'test_GD_1atm_3D'
file = gd.de.Make_dir(save_path)
to_optimize = 'energy_demand'   # key in data from de.get_data()
searches = 12    # number os total searches
trial_points = 240  # number of trial start_points. best ones will be used for searches
~~~

The save path for the generated CSV files is also set here, under `save_path`. The example will first evaluate `trial_points` number of random points in a small scale bruteforce study. The best `searches` results will be used as the starting points for gradient descent based optimization. 

Then in the dictionary called `ranges`, you can set a range for all control parameters (except gases and fractions), souch as `R_E = [0.5e-6, 20e-6]`. Now the min value for R_E is 0.5 um, and the maximal value will be 20 um. If you don't want to change a control parameter, you can set a fixed value as `R_E = [5.0e-6]`. Now R_E will be 5 um, and the algorithm won't optimize it.

Remember, that a search will converge to the nearest local minimum, that's why you need multiple. A search may include several simulations, so you may need to wait for the first results to be displayed. The searches are independent, so you can sut down the program, before it finishes. Afterward, some plots help you determine, if the searches were convergent.

### Gradient_descent_2D.ipyb

This example can read CSV files containing a 2D bruteforce parameter sweep, and display the energy demand on a contour plot like hte one above. Linear interpolation is used to fill missing values. The example can also visualize a 2D gradient descent based search.

<a name="create_plots"/>

## Create plots

Example file: [Create_plots.ipynb](Create_plots.ipynb)

|![Image not found](https://user-images.githubusercontent.com/42745647/215842722-d057eefe-c54f-4046-ab8d-02a80ab171ef.png)|
|:---:|
|A plot showcasing the energy demand of ammonia production as a function of equilibrium bubble size and ambient pressure.|

This example shows you how to create simple plots with [matplotlib.pyplot](https://matplotlib.org/), like the one above. In the plots, the energy demand is examined with respect to one ore more control parameters.
