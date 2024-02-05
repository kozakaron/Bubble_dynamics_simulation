
# Bubble dynamics simulation documentation

This document is automatically generated from in code documentation. (multiline comments under function headers) 

Date of generation: 2024.02.05 20:44:32 (YYYY.MM.DD HH:MM:SS) 
## Table of contents

0. [**inp_data_extractor.py**](#bookmark_inp_data_extractor)
1. [**excitation.py**](#bookmark_excitation)
2. [**full_bubble_model.py**](#bookmark_full_bubble_model)
3. [**gradient_descent.py**](#bookmark_gradient_descent)
4. [**data.py**](#bookmark_data)
5. [**doc_generator.py**](#bookmark_doc_generator)
<a name="bookmark_inp_data_extractor"/>

## inp_data_extractor.py

[Link to file](inp_data_extractor.py)

~~~
This program extracts data from .inp files and creates parameters.py
Recommended usage:
    importing: import inp_data_extractor as inp
    usage: inp.extract(path)
~~~

Internal dependencies: 
 * data.py

External dependencies: 
 * importlib
 * numpy
 * os

### Functions

* **extract**

	~~~Python
	def extract(path)
	~~~

	~~~
	Extracts data from .inp file and creates parameters.py in the current working directory. Arguments:
	 * path (str): path to .inp file
	~~~

* **print_array**

	~~~Python
	def print_array(array, width=0, comments=[], columns=[], max_len=0)
	~~~

	~~~
	Prints a 1D or 2D array or list. Can leave comments above columns and after lines. Arguments:
	 * array (list or np.ndarray): array to be printed
	 * width (int): width of each element in characters
	 * comments (list): list of comments after each line
	 * columns (list): list of comments above each column
	 * max_len (int): maximum number of elements in one line, if 0, all elements are printed in one line
	~~~


### Global variables

Actual values might differ
 * **comment**: *str* = `'!'`
<a name="bookmark_excitation"/>

## excitation.py

[Link to file](excitation.py)

~~~
This file contains different pressure excitation functions. 
Specify which one to use under "___Settings___" in full_bubble_model.py 

Use getExcitation(excitation_type) to get the excitation function and its arguments. 
Use plotExcitation(excitation_type) to plot the excitation function. 
~~~


External dependencies: 
 * numpy
 * os
 * matplotlib.pyplot

### Functions

* **getExcitation**

	~~~Python
	def getExcitation(excitation_type='no_excitation')
	~~~

	~~~
	Returns the pressure excitation function and the arguments it takes. Use plotExcitation() to plot the excitation function.
	Available excitation types:
	 * 'no_excitation': constant ambient pressure
	 * 'two_sinusoids': two sinusoids with different frequencies and amplitudes, and a phase shift between them
	 * 'sin_sqr': sinusoid squared with only n amplitude cycles
	 * 'slow_expansion': decrease from abient pressure to min_pressure (decay_time), followed by a slow increase back to ambient pressure (increase_time)
	 * 'sin_impulse_flat_ends': sinusoid with only n amplitude cycles, the ends are smoothed out
	 * 'sin_impulse': sinusoid with only n amplitude cycles, the ends are not smoothed out
	 * 'sin_impulse_logf': same as 'sin_impulse', but log10(freq) is used instead of freq
	 * 'double_sin_impulse': n cycle of two sinusoids with different frequencies and amplitudes and no phase shift; the 2nd frequency is freq_ratio times the 1st frequency
	 * 'multi_sin_impulse': 5 cycles of sinusoids with different frequencies and amplitudes for each cycle
	 * 'double_multi_sin_impulse': 5 cycles of two sinusoids with different frequencies and amplitudes for each cycle; 2 sine waves are used each time, similar to 'double_sin_impulse'
	 
	Returns:
	 * Excitation: jitted excitation function
	    	* inputs: t (time, sec), P_amb (ambient pressure, Pa), args (list of control parameters)
	    	* outputs: p_Inf (excitation pressure, Pa), p_Inf_dot (excitation pressure derivative, Pa/s)
	 * args: list of argument names the Excitation function takes
	 * units: list of units of the arguments
	 * defaults: list of default values for the arguments
	~~~

* **plotExcitation**

	~~~Python
	def plotExcitation(excitation_type, ex_args, t_int=[0.0, 0.0001], base_name='', format='png', presentation_mode=False)
	~~~

	~~~
	Plots the excitation function.
	Arguments:
	 * excitation_type: string, excitation type
	 * ex_args: list of control parameters
	 * t_int: time interval to plot excitation (default: [0, 0.0001] [s])
	 * base_name: save plots as image (default: '' alias do not save) | 
	           use base_name='plot' --> plot_1.png, plot_2.png |  
	           use base_name='images/plot' to save into images folder |  
	           using a folder for images is recommend |  
	           this folder have to be created manually
	 * format: format of the saved images (available: png, pdf, ps, eps, svg)
	 * presentation_mode: if True, the plot will be in presentation mode (default: False)
	~~~

<a name="bookmark_full_bubble_model"/>

## full_bubble_model.py

[Link to file](full_bubble_model.py)

~~~
This program contains the JIT compiled (high performance) ODE function (_f()), the numerical solvers (solve()), and several supplementary functions. 
Before use, make sure to set the variables below "___Settings___" to your needs. 

Usage:
 * Use inp_data_extractor.py to turn a .inp file into parameters.py
 * Set the parameters in the "___Settings___" section
 * Import this file: from Bubble_dynamics_simulation import full_bubble_model as de
 * Assemble the control parameter dictionary: cpar = de.example_cpar() | 
   You may print cpar: de.print_cpar(cpar) | 
   The result can be copied, modified and used as code. 
 * You may solve the differential equation: num_sol, error_code, elapsed_time = de.solve(cpar)
 * Retrieve post processing data: data = de.get_data(cpar, num_sol, error_code, elapsed_time)
 * Use de.plot(cpar) to plot the results. This function has several costumization options and extra plots.

 Notes:
  * See the documentation or run help(de.solve) for more information. Substitute any function name instead of de.solve.
  * Use the Make_dir class to save simulation results.
  * See the example files. Example files start with a capital letter, and have the format .ipynb
  * You may define any excitation function in excitation.py with a fix number of excitation control parameters.
  * You can plot any extra variables with plot() by setting plot_extra=True. 
    To do this, modify the bottom lines in the _f() function, and the first lines of plot().
~~~

Internal dependencies: 
 * excitation.py
 * parameters.py

External dependencies: 
 * importlib
 * numpy
 * os
 * matplotlib.pyplot
 * psutil
 * socket
 * time

### Classes

* **Make_dir**

	~~~Python
	class Make_dir(self, folder_name, file_base_name='output_', separator=',')
	~~~

	~~~
	Class for saving simulation data to csv files. Constructor arguments:
	 * folder_name: name of the folder to save the csv files into (e.g. 'folder', 'folder/subfolder')
	 * file_base_name: base name of the csv files (default: 'output_' --> 'output_1.csv', 'output_2.csv', ...)
	 * separator: separator character in the csv file (default: ',')
	~~~

	* *close*

		~~~Python
		def close(self)
		~~~

		~~~
		Closes the currently opened csv file.
		~~~

	* *new_file*

		~~~Python
		def new_file(self)
		~~~

		~~~
		Creates a new csv file from file_base_name and an unique number. (e.g. 'output_1.csv', 'output_2.csv', ...)
		Automatically closes opened file. Use write_line() to write data into the file.
		~~~

	* *write_line*

		~~~Python
		def write_line(self, data)
		~~~

		~~~
		Writes the data dict into the currently opened csv as a new line. The line contains the values of the keys in the data dict.
		~~~

	* *write_solution*

		~~~Python
		def write_solution(self, data, num_sol, file_base_name)
		~~~

		~~~
		Writes the data dict and the numerical solution into two new csv files. Arguments:
		 * data: dictionary containing the simulation data from get_data()
		 * num_sol: solution of the differential equation from solve()
		 * file_base_name: base name of the csv files (e.g. 'name' --> 'name_data.csv', 'name_num_sol.csv')
		~~~

	* *write_string*

		~~~Python
		def write_string(self, string, file_base_name)
		~~~

		~~~
		Writes the string into a new txt file. Also saves header with get_settings_and_info(). Arguments:
		 * string: arbitrary string to write into the file
		 * file_base_name: base name of the txt file (e.g. 'name' --> 'name.txt')
		~~~

* **dotdict**

	~~~Python
	class dotdict()
	~~~

	~~~
	Dot notation access to dictionary attributes. 
	Instead of dictionary['key'] you can use dictionary.key
	~~~


### Functions

* **VapourPressure**

	~~~Python
	def VapourPressure(T)
	~~~

	~~~
	Legacy function, use vapour_pressure(T) instead
	~~~

* **Viscosity**

	~~~Python
	def Viscosity(T)
	~~~

	~~~
	Legacy function, use viscosity(T) instead
	~~~

* **copy**

	~~~Python
	def copy(input)
	~~~

	~~~
	Deep copies the input. Use: copied = copy(original). 
	Deep copy means that the input and the output are independent from each other, nothing is copied by reference.
	~~~

* **example_cpar**

	~~~Python
	def example_cpar(normal_dict=False)
	~~~

	~~~
	Provides an example of the control parameter dictionary. Use print_cpar() to print it. Parameters:
	 * normal_dict: if True, returns a normal dictionary, else returns a dotdict
	 
	 Returns:
	 * cpar: control parameter dictionary
	~~~

* **get_data**

	~~~Python
	def get_data(cpar, num_sol, error_code, elapsed_time)
	~~~

	~~~
	This function gets the numerical solution and the control parameters, and returns some datas about the simulation. Arguments:
	 * cpar: control parameters (dict or dotdict)
	 * num_sol: numerical solution from solve()
	 * error_code: error code from solve()
	 * elapsed_time: elapsed time from solve()
	
	Returns:
	 * data: dotdict with the post processing data (e.g. collapse time, energy demand, etc.)
	~~~

* **get_errors**

	~~~Python
	def get_errors(error_code, printit=False)
	~~~

	~~~
	 * Input: error_code (int) 
	 * Output: list of error codes (str)
	 * Also prints colored errors, if printit=True
	~~~

* **get_settings_and_info**

	~~~Python
	def get_settings_and_info()
	~~~

	~~~
	Returns a string with information about the computer, and the settings of the simulation.
	~~~

* **plot**

	~~~Python
	def plot(cpar, t_int=[0. 1.], n=5.0, base_name='', format='png', LSODA_timeout=30.0, Radau_timeout=300.0, presentation_mode=False, plot_pressure=False, plot_extra=False, show_legend=False, show_cpar=True)
	~~~

	~~~
	This funfction runs solve(), get_data(), print_data(), and plots the numerical solution. 
	By default, R(t) and T(t) are plotted on the first plot. 
	The amount of chemical species are plotted on the second plot. 
	Optionally, the internal and external pressure can be plotted on the third plot. 
	Parameters:
	 * cpar: control parameters in a dictionary
	 * t_int: time interval to solve the diffeq in (default: [0, 1] [s]) |   
	       graphs will be plotted in this intervall, if not default
	 * n: how long should the plotted time interval be compared to the collapse time (default: 5 [-])
	 * base_name: save plots as image (default: '' alias do not save) | 
	           use base_name='plot' --> plot_1.png, plot_2.png |  
	           use base_name='images/plot' to save into images folder |  
	           using a folder for images is recommend |  
	           this folder have to be created manually
	 * format: format of the saved images (available: png, pdf, ps, eps, svg)
	 * LSODA_timeout, Radau_timeout: timeout (maximum runtime) for different solvers in solve() in seconds
	 * presentation_mode: if True, the plot will be in presentation mode (default: False)
	 * plot_pressure: if True, the internal and external pressure will be plotted (default: False)
	 * plot_extra: if True, extra dimensions will be plotted (default: False) | 
	               To use this, you have to modify the end of the _f() function and set extra_dims and extra_dim_labels below. 
	 * show_legend: if True, the legend will be visible with every single species (default: False)
	 * show_cpar: if True, the control parameters will be printed on the plot (default: False)
	~~~

* **print_cpar**

	~~~Python
	def print_cpar(cpar, without_code=False, print_it=True)
	~~~

	~~~
	Prints the control parameters (cpar) in an organised way. Arguments:
	 * cpar: control parameters (dict or dotdict)
	 * without_code: if True, an easier to read version is printed. If False, then the result is a valid python code
	 * print_it: if True, the function will print the text. If False, it will return with the text (sting)
	~~~

* **print_data**

	~~~Python
	def print_data(cpar, data, print_it=True)
	~~~

	~~~
	Prints the data dictionary from get_data() in an organised way. Arguments:
	 * cpar: control parameters
	 * data: data dictionary from get_data()
	 * print_it: if True, the function will print the text. If False, it will return with the text (sting)
	~~~

* **simulate**

	~~~Python
	def simulate(kwargs)
	~~~

	~~~
	This function runs solve() and get_data(), then return with data. 
	Input and output is (or can be) normal dictionary. 
	It is used for multithreading (e.g. in bruteforce_parameter_study.inp). 
	The input (kwargs) is a dictionary with the keyword-argument pairs of solve().  
	~~~

* **solve**

	~~~Python
	def solve(cpar, t_int=[0. 1.], LSODA_timeout=30.0, Radau_timeout=300.0, extra_dims=0)
	~~~

	~~~
	This funfction solves the differential equation, and returns the numerical solution.
	Parameters:
	 * cpar: control parameters
	 * t_int: time interval
	 * LSODA_timeout: timeout for LSODA solver in seconds
	 * Radau_timeout: timeout for Radau solver in seconds
	 * extra_dims: add extra dimensions to the initial condition array (initial value: 0.0) | 
	               use it to plot extra variables (e.g. energy) during the simulation
	
	Returns:
	 * num_sol: numerical solution. Use num_sol.t and num_sol.y to get the time and the solution. Can be None
	 * error_code: see de.error_codes: dict, de.get_errors()
	 * elapsed_time: elapsed time
	~~~

* **vapour_pressure**

	~~~Python
	def vapour_pressure(T)
	~~~

	~~~
	Calculates the vapour pressure of water as a function of temperature (T) in Kelvin. 
	~~~

* **viscosity**

	~~~Python
	def viscosity(T)
	~~~

	~~~
	Calculates the dynamic viscosity of water as a function of temperature (T) in Kelvin. Pressure dependence is neglected. 
	~~~


### Global variables

Actual values might differ
 * **enable_dissipated_energy**: *bool* = `True`
 * **enable_evaporation**: *bool* = `False`
 * **enable_heat_transfer**: *bool* = `True`
 * **enable_reactions**: *bool* = `True`
 * **error_codes**: *dict* = `{'xx0': {'describtion': 'succecfully solved with LSODA solver', 'color': 'green'...`
 * **excitation_args**: *list* = `['p_A', 'freq', 'n']`
 * **excitation_defaults**: *list* = `[-200000.0, 30000.0, 1.0]`
 * **excitation_type**: *str* = `'sin_impulse'`
 * **excitation_units**: *list* = `['Pa', 'Hz', '-']`
 * **keys**: *list* = `['ID', 'R_E', 'ratio', 'P_amb', 'alfa_M', 'T_inf', 'P_v', 'mu_L', 'rho_L', 'gase...`
 * **target_specie**: *str* = `'NH3'`
<a name="bookmark_gradient_descent"/>

## gradient_descent.py

[Link to file](gradient_descent.py)

~~~
This module is used to optimize the energy demand by automatically adjusting the control parameters. 
A search (using search() or gradient_descent()) can start in any point (e.g. in a random point from rand_point()). 
In each step, the gradient is calculated with finite differences. The step is opposite to the normed gradient. 
The step size is automatically controlled. To set ranges and arguments, see documentation_for_gradient_descent.pdf
~~~

Internal dependencies: 
 * full_bubble_model.py

External dependencies: 
 * importlib
 * numpy
 * random
 * time

### Functions

* **evaluate**

	~~~Python
	def evaluate(point, to_optimize, t_int, LSODA_timeout, Radau_timeout)
	~~~

	~~~
	Runs the simulation, and returns with data extended with output (what we want to optimize). Arguments:
	 * point: dict, control parameter which we want to evaluate
	 * to_optimize: str, name of the output we want to optimize (e.g. 'energy_demand')
	 * t_int, LSODA_timeout, Radau_timeout: see de.solve() for more info
	 
	 Returns:
	  * data: dict, simulation results from de.get_data()
	  * success: bool, True if the simulation was successful, False otherwise
	~~~

* **evaluate_kwargs**

	~~~Python
	def evaluate_kwargs(kwargs)
	~~~

	~~~
	Call evaluate() with a dictionary containing the arguments. Arguments:
	 * kwargs: dict, containing the arguments for evaluate()
	     
	 Returns:
	 * data: dict, simulation results from de.get_data()
	 * point: dict, control parameter which we want to evaluate
	~~~

* **gradient_descent**

	~~~Python
	def gradient_descent(ranges, path, to_optimize, start_point, step_limit=100, first_step=0.01, min_step=0.0001, delta=1e-06, verbose=True, t_int=[0. 1.], LSODA_timeout=30, Radau_timeout=300)
	~~~

	~~~
	Gradient search. Starts at start_point, always go in the direction of the local gradient. Step size decays, if the output at the next step would
	bigger then at the current step, or if too many steps were taken without decay (to avoid back &forth stepping). Search ends, if the step_size
	decays to be smaller than min_step*interval_width, or if gradient fails repeatedly.     Arguments:
	 * ranges: dict, ranges of the parameters ([single_value] or [min, max])
	 * path: str, save location. If None or '', datas will be returned (not recommended due to memory usage), otherwise last data will be returned
	 * to_optimize: str, name of the output we want to optimize (e.g. 'energy_demand')
	 * start_point: dict, cpar where the search starts
	 * step_limit: int, maximum number of steps
	 * first_step: float, first step size
	 * min_step: float, minimum step size
	 * delta: float, step size for the finite difference
	 * t_int, LSODA_timeout, Radau_timeout: see de.solve() for more info
	 * verbose: bool, print stuff
	
	Returns:
	 * datas: list of dicts, simulation results from de.get_data() if path is None or last data if path is given
	 * last_bests: list of floats, list of the best outputs in each step
	 * elapsed_time: float, elapsed time in seconds
	~~~

* **rand_point**

	~~~Python
	def rand_point(ranges, ID=0, padding=0.001)
	~~~

	~~~
	Generates a random cpar dict within ranges. Arguments:
	 * ranges: dict, ranges of the parameters ([single_value] or [min, max])
	 * ID: int, ID for the point
	 * padding: float, padding for the ranges (0.0 <= padding <= 0.5) to avoid points too close to the edge
	
	Returns:
	 * point: dict, cpar with random values (warning: not a dotdict)
	~~~

* **search**

	~~~Python
	def search(kwargs)
	~~~

	~~~
	Call gradient_descent() with a dictionary containing the arguments.
	~~~

<a name="bookmark_data"/>

## data.py

[Link to file](data.py)

~~~
This python file contains numerical data missing from the .inp OpenSmoke files.
Recommended usage:
   importing: import data
   usage: print('Molar mass of hydrogen: =', data.W['H'])
~~~


### Functions

* **calculate_missing_constants**

	~~~Python
	def calculate_missing_constants()
	~~~


### Global variables

Actual values might differ
 * **W**: *dict* = `{'H': 1.00797, 'HE': 4.0026, 'LI': 6.941, 'BE': 9.01218, 'B': 10.81, 'C': 12.011...`
 * **header**: *str* = `'""", This is an automatically generated file., This python file contains all the n...'`
 * **lambdas**: *dict* = `{'HE': 0.151, 'NE': 0.0491, 'AR': 0.0177, 'KR': 0.00943, 'H2': 0.1805, 'O2': 0.0...`
 * **physical_constants**: *dict* = `{'c_L': {'value': 1483.0, 'comment': 'Liquid sound speed at 30 °C [m/s]'}, 'rho_...`
 * **valid_elements**: *str* = `'H, HE, LI, BE, B, C, N, O, F, NE, NA, MG, AL, SI, P, S, CL, AR, K, CA, SC, TI, V...'`
<a name="bookmark_doc_generator"/>

## doc_generator.py

[Link to file](doc_generator.py)

~~~
This software is used to automatically generate the documentation for the project. Function declarations and in code documentation is used for this porpuse. 
You should doucument a function directly underneath the function header in a multiline comment:

def function_name(arg1, arg2):
"""Comment goes here. Define the arguments and the return valuea as well."""
    pass

Use: just run the file. A markdown called DOCUMENTATION.md will be automatically generated. 
You may specify, which modules to document by adding them to the to_document list. 
Set document_private to True to include private functions and variables in the documentation. (e.g. _f())
~~~

Internal dependencies: 
 * data.py
 * full_bubble_model.py
 * doc_generator.py
 * excitation.py
 * gradient_descent.py
 * inp_data_extractor.py

External dependencies: 
 * os

### Functions

* **list_all**

	~~~Python
	def list_all(module, show_private=False)
	~~~

	~~~
	Returns a list of all modules, classes, functions and variables in the given module. 
	Can hide private elements. (starting with '_') 
	Submodules are marked with a '#' in front of their name. 
	~~~

* **main**

	~~~Python
	def main()
	~~~

	~~~
	Generates the documentation for the project.
	~~~

* **print_function**

	~~~Python
	def print_function(fun, indent=0)
	~~~

	~~~
	Returns the function declaration as a string. 
	The string is formatted to be a markdown code block, and indented by the given number of tabs.
	~~~

* **print_help**

	~~~Python
	def print_help(fun, indent=0)
	~~~

	~~~
	Returns the __doc__ of the function as a string. 
	The string is formatted to be a markdown code block, and indented by the given number of tabs.
	~~~

* **print_variable**

	~~~Python
	def print_variable(name, value, use_type=False)
	~~~

	~~~
	Returns the name and value of a variable as a formatted string.
	~~~


### Global variables

Actual values might differ
 * **document_private**: *bool* = `False`
 * **header**: *str* = `', # Bubble dynamics simulation documentation, , This document is automatically gene...'`
 * **to_document**: *list* = `[<module 'inp_data_extractor' from 'd:\\parameter_studies\\Bubble_dynamics_simul...`
