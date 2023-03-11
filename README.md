#  Bubble_dynamics_simulation

Tutorial video: [Magyar](https://youtu.be/YWsT1ktUzVw), [English]()

Author: Kozák Áron

Budapest University of Technology and Economics, 
Department of Hydrodynamic Systems

email: kozi0223@gmail.com

## Table of content

1. [**Python basics**](#python_basics): Feel free to skip this part. Helps with installations, briefly introduces some special python construct, and links tutorials. Might be helpful, if you are new to python.

2. [**The simulation**](#the_simulation): Introduces some of the basic functions in *full_bubble_model.py*. Shows, how to solve the differential equation with different control parameters. Shows how to plot, process, or save the results.

3. [**Bruteforce parameter sweep**](#bruteforce_parameter_sweep): A simple multithread example. A control parameter space is sweeped with the given resolution. Results are saved into CSV files.

4. [**Read CSV files**](#read_csv_files): This example showcases how to read CSV files from a given folder. The data is loaded into a pandas dataframe. Different statistics and data manipulation can be done with the dataframe.

5. [**Create plots**](#create_plots): This example is about creating simple 1D plots using the simulation.

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
* **.py** files: The original python format. It's just plain text, like *.c* *.cpp* *.h*. You can edit them in any text editor, and run them from a command prompt without any IDEs. Only *.py* files can be imported into other projects as libraries, this is the reason I use them.
* **.ipynb** files: The notebook format, in a text editor you will see a html like structure, which has to be rendered, thus you need to use an IDE. A notebook contains code cells, which can run separately, and have they own output. There are also markdown cells, which displays formatted text with headers, latex formulas, images, etcetera. Some IDEs can generate a table of content from the headers.

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

Also there is a shortcut to create lists:
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

Example file: *full_bubble_mode.ipynb* <br>

For import: *full_bubble_mode.py*, *parameters.py* <br>

The mathematical model is detailed in *full_bubble_mode.ipynb*. This is a notebook, where each function has a little demo, and all the formulas are attached in markdown cells. In most editors, you can navigate the file with the markdown headers. Note, that notebooks can not be imported in other codes, therefore all the code cells are copied to *full_bubble_mode.py*. If you modify the notebook (.ipynb), make sure to copy the changes into the python file (.py).

###  Parameters, automatic generation

For import: *_inp_data_extractor.py*, *data.py*, a *.inp* file <br>

All the numerical constants and coefficients are stored in a separate file. This is an automatically generated file, always named *parameters.py*, and it's essential for the simulation. You can import this file and use it's variables:

~~~Python
import parameters as par
par.model # get the name of the .inp file, that was used to create parameters.py
~~~

~~~
chemkin_AR_HE
~~~

Make sure, you don't change the variables in par in the code! If you `import full_bubble_model as de`, you can acces the parameters, as `de.par` . <br>
You can automatically generate *parameters.py* from a .inp file, such as one from the *INP file examples* folder. To do this, you will need *inp_data_extractor.py* and *data.py*. In *data.py* you will find all constants missing from the .inp files, like molar masses and lambda values. For the generation, run the following code:

~~~Python
# Directory for .inp file:
path = 'INP file examples\\chemkin_AR_HE.inp'   # always use \\ instead of \

import inp_data_extractor as inp
inp.extract(path)   # this command creates parameters.py, it owerwrites the existing one
~~~

The `inp.extract()` function might give you colored messages. The blue ones can be neglected, the yellow ones should be considered, but won't result in an error, but the red ones will cause errors, and must be investigated.

Now you can `import parameters as par`. Note, that if you, or your code overwrites a file, you have to restart the kernel to see the changes. Alternatively, you can use `import importlib`, then you `import parameters as par`, and when you change the *parameters.py* file, you just have to run `importlib.reload(par)` to use the changes.

###  Control parameters

The control parameters are stored in a dictionary for easier managing. (Not in the high performance part) An example of such a dictionary:

~~~Python
import full_bubble_model as de

cpar = de.dotdict(dict(
	ID=0,
	R_E=1e-6, # [m]
	ratio=5.0, # [-]
	P_inf=25e5, # [Pa]
	alfa_M=0.05, # [-]
	T_inf=303.15, # [K]
	surfactant=1.0, # [-]
	gases=[de.par.index['AR']], # indexes of species in initial bubble
    	fractions=[1.0], # molar fractions of species in initial bubble
))
cpar.P_v = de.VapourPressure(cpar.T_inf) # [Pa]
cpar.mu_L = de.Viscosity(cpar.T_inf) # [Pa*s]
~~~

This cpar will be passed to most of the not JIT-ted functions. Don't forget, that cpar most contain an ID number, and all 10 control parameters, including the saturated vapour pressure and the viscosity of the liquid, which is calculated from the temperature, otherwise not independent. You can of course use a constant value for them. Note, that viscosity is also pressure dependent, which is neglected here, and this function only gives the viscosity of water.

###  Simulating

You should start with

~~~Python
import full_bubble_model as de
~~~

The most important functions: <br>

* **plot(cpar, *t_int*, *n*, *base_name*, *LSODA_timeot*, *Radau_timeout*)**: This funfction solves the differential equation, and plots it <br>

	Arguments:

	* cpar

	* t_int (optional): time interval, the default is [0, 1.0] seconds. Graphs will be plotted in this intervall, if you change it.

	* n (optional): how long should the plotted time interval be compared to the collapse time, default: 5 [-]

	* base_name (optional): save plots as .png, default: don't save. Use `base_name='plot'` --> plot_1.png, plot_2.png. Use `base_name='images\\plot'` to save into images folder. Using a folder for images is recommended. This folder have to be created manually

	* LSODA_timeout, Radau_timeout (optional): maximum time the different solvers are allowed to run. After timeout seconds, the solver will terminate, and return with an error code. The default is 30 [s] for LSODA and 300 [s] for Radau.

	Returns: - <br>

	Example:

	~~~Python
	de.plot(cpar)
	~~~

	~~~
	succecfully solved with LSODA solver
	~~~

![image](https://user-images.githubusercontent.com/42745647/215813926-de63814e-8c04-4dbb-8e00-c65ee00e135b.png)
![image](https://user-images.githubusercontent.com/42745647/215814012-e4f4cdd4-29c1-4420-8c0a-07d6bb9093fb.png)

~~~
Control parameters:
    ID = 11
    R_E = 1.00 [um]
    ratio = 5.00 [-]
    P_inf = 25.00 [bar]
    alfa_M = 0.05 [-]
    T_inf = 30.00 [°C]
    P_v = 4245.13 [Pa]
    mu_L = 0.81 [mPa*s]
    surfactant = 1.00 [-]
    100% AR, 0% HE
Simulation info:
    error_code = 0
    elapsed_time = 1.10 [s]
    steps = 9539 [-]
Final state:
    R_final = 1.01 [um];   R_dot_final =1.129169540508806e-17 [m/s];   T_final = 303.15 [K]
    n_H2 =5.826818607281378e-17 [mol]; n_O2 =2.864245341970364e-17 [mol]
    Final molar concentrations: [mol/cm^3]
        H: -4.994195859260988e-21;  H2: 1.3632275940074182e-05;  O: -3.4107975464517105e-15;  O2: 6.701115221438125e-06;  
        OH: -1.0036204724547934e-19;  H2O: 1.6842204457612123e-06;  N2: 0.0;  HO2: -2.4113237248546457e-11;  
        H2O2: 2.300816704644442e-07;  AR: 0.0010263315529529951;  HE: 0.0;  OHEX: 8.056090647320269e-45;  
        
Results:
    collapse_time = 9.422505710336362e-08 [s]
    T_max = 4048.15 [K]
    expansion work = 1.264638360255316e-09 [J]
    hydrogen production = 10766.07 [MJ/kg]
~~~

* **solve(cpar, *t_int, LSODA_timeout, Radau_timeout*)**: This funfction solves the differential equation, and returns the numerical solution. Use, if you want to get the numerical solution itself for further investigations. <br>

	Arguments:

	* cpar

	* t_int (optional): time interval, the default is [0, 1.0] seconds.

	* LSODA_timeout, Radau_timeout (optional): maximum time the different solvers are allowed to run. After timeout seconds, the solver will terminate, and return with an error code. The default is 30 [s] for LSODA and 300 [s] for Radau.

	Returns:

	* num_sol: an object containing the numerical solutions. Time is stored in array num_sol.t, the solution in array num_sol.y

	* error_code: a number from 0 to 6. If 4 or greater, the solvers weren't succesful.

	* elapsed_time: total runtime of the function in seconds.

	Example:

	~~~Python
	num_sol, error_code, elapsed_time = de.solve(cpar)
	~~~

* **get_data(cpar, num_sol, error_code, elapsed_time)**: This function gets the numerical solution and the control parameters, and returns some datas about the simulation. <br>

	Arguments:

	* cpar

	* num_sol: from solve()

	* error_code: from solve()

	* elapsed_time: from solve()

	Returns: A dictionary containing several inforamtions about the simulation, souch as the control parameters, error code, final results.

	Example:

	~~~Python
	data = de.get_data(cpar, num_sol, error_code, elapsed_time)
	data
	~~~

	~~~
	{'ID': 0,
	 'R_E': 3e-06,
	 'ratio': 20.0,
	 'P_inf': 5000000.0,
	 'alfa_M': 0.05,
	 'T_inf': 303.15,
	 'P_v': 4245.12571625229,
	 'mu_L': 0.0008148611589373901,
	 'surfactant': 0.25,
	 'gases': [3],
	 'fractions': [1.0],
	 'error_code': 0,
	 'elapsed_time': 2.483180046081543,
	 'steps': 17628,
	 'collapse_time': 8.006066235613211e-07,
	 'T_max': 3218.745891558094,
	 'x_final': array([ 3.00142939e-06,  3.25564626e-17,  3.03150000e+02,  1.05099081e-33,
			            8.78391179e-08,  4.66028275e-31,  1.98129278e-03,  2.46247464e-31,
			            1.68422045e-06,  0.00000000e+00,  5.89285538e-12,  5.40273169e-06,
			            0.00000000e+00,  0.00000000e+00, -2.25892096e-71]),
	 'n_H2': 9.9485770654553e-18,
	 'n_O2': 2.2439938309400862e-13,
	 'work': 4.515209424703302e-06,
	 'energy': 225133087.80139184}
	~~~

* **print_data(data, *print_it*)**: This function prints the data dictionary in an organised way. <br>

	Arguments:

	* data: from get_data()

	* print_it (optional): Default is True, in this case, the function prints the text to the console. If set to False, the function will return the text, and won't print anything.

	Returns: - <br>

	Example:

	~~~Python
	de.print_data(data)
	~~~

	~~~
	Control parameters:
	    ID = 0
	    R_E = 3.00 [um]
	    ratio = 20.00 [-]
	    P_inf = 50.00 [bar]
	    alfa_M = 0.05 [-]
	    T_inf = 30.00 [°C]
	    P_v = 4245.13 [Pa]
	    mu_L = 0.81 [mPa*s]
	    surfactant = 0.25 [-]
	    100% O2
	Simulation info:
	    error_code = 0
	    elapsed_time = 2.48 [s]
	    steps = 17628 [-]
	Final state:
	    R_final = 3.00 [um];   R_dot_final =3.255646260115485e-17 [m/s];   T_final = 303.15 [K]
	    n_H2 =9.9485770654553e-18 [mol]; n_O2 =2.2439938309400862e-13 [mol]
	    Final molar concentrations: [mol/cm^3]
	        H: 1.0509908121293901e-33;  H2: 8.783911792457889e-08;  O: 4.660282747441842e-31;  O2: 0.0019812927762544595;  
	        OH: 2.4624746417526286e-31;  H2O: 1.6842204457221403e-06;  N2: 0.0;  HO2: 5.892855382810517e-12;  
	        H2O2: 5.402731685915859e-06;  AR: 0.0;  HE: 0.0;  OHEX: -2.2589209631942755e-71;  
	        
	Results:
	    collapse_time = 8.006066235613211e-07 [s]
	    T_max = 3218.75 [K]
	    expansion work = 4.515209424703302e-06 [J]
	    hydrogen production = 225133087.80 [MJ/kg]
	~~~

* **simulate(kwargs)**: This function runs solve() and get_data(), then returns with data. Both the input and the output is (or can be) a normal dictionary. This function is used for multithreading in the bruteforce parameter sweep and in other examples. <br>

	Arguments:

	* kwargs: a dictionary containing all the keyword arguments, you would pass to `solve()`.

	Returns:

	* data: same, as in get_data(), but normal dictionary

	Example:

	~~~Python
	data = de.simulate( dict(cpar=cpar, t_int=np.array([0.0, 1.0])) )
	~~~

* **Make_dir**: class for saving things into CSV files. To use, first run:

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

	You can close the CSV to save it:

	~~~Python
	file.close()
	~~~

	You can also save the numerical solution:

	~~~Python
	num_sol, error_code, elapsed_time = de.solve(cpar)
	data = de.get_data(cpar, num_sol, error_code, elapsed_time)
	file.write_solution(data, num_sol, 'testfile')
	~~~

	Two files will be created in the *test* folder. *testfile_data.csv* contains the data dictionary, while *testfile_num_sol.csv* contains the numerical solution. This function is independent of new_file() and close().

<a name="bruteforce_parameter_sweep"/>

##  Bruteforce parameter sweep

Example file: *bruteforce parameter sweep.ipynb* <br>

For import: *full_bubble_mode.py*, *inp_data_extractor.py*, *data.py*, a *.inp* file <br>

In this example, you can create a parameter study. You give a list for all control parameters, containing all their possible values. Then the program makes a long list of all the possible combinations. Affterward importing [multiprocessing](https://docs.python.org/3/library/multiprocessing.html#using-a-pool-of-workers), you can solve the simulation for all combinations using all your computer cores. (e.g. R_E=[1, 2, 3, 4] [um]; expansion ratio=[5, 10, 15, 20] [-] --> total 4*4=16 combinations) <br>

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

Normally, using a for loop, it would take 16 seconds to run all the combinations, but using all your cpu threads, you can get it done much faster. The example might seem a bit more complicated, but the core mechanic is the same. The examlpe uses pool.imap_unordered, instead of pool.imap to get result faster, but in arbitrary order. The example uses the Make_dir class to save the results to .csv files.

<a name="read_csv_files"/>

##  Read CSV files

Example file: *read csv files.ipynb* <br>

For import: *full_bubble_mode.py*, *inp_data_extractor.py*, *data.py*, a *.inp* file <br>

After running a parameter study, you will have several .csv files with the results. This example uses [pandas](https://pandas.pydata.org/docs/), a popular data science library. The program lists and loads all csv files from the given folder into a dataframe, and makes some data manipulations, like filtering solutions, and also gets some statistics. <br>

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

<a name="create_plots"/>

## Create plots

Example file: *create plots.ipynb* <br>

For import: *full_bubble_mode.py*, *inp_data_extractor.py*, *data.py*, a *.inp* file <br>

![image](https://user-images.githubusercontent.com/42745647/215842722-d057eefe-c54f-4046-ab8d-02a80ab171ef.png)

This example shows you how to create simple plots with matplotlib.pyplot, like the one above.
