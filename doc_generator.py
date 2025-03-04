"""
This software is used to automatically generate the documentation for the project. Function declarations and in code documentation is used for this porpuse. 
You should doucument a function directly underneath the function header in a multiline comment:

def function_name(arg1, arg2):
\"\"\"Comment goes here. Define the arguments and the return valuea as well.\"\"\"
    pass

Use: just run the file. A markdown called DOCUMENTATION.md will be automatically generated. 
You may specify, which modules to document by adding them to the to_document list. 
Set document_private to True to include private functions and variables in the documentation. (e.g. _f())
"""

# To use:
from datetime import datetime
import os

# To document
import inp_data_extractor as inp
import excitation
import full_bubble_model as de
import gradient_descent as gd
import data
import doc_generator

to_document = [inp, excitation, de, gd, data, doc_generator]
document_private = False

header = f'''
# Bubble dynamics simulation documentation

This document is automatically generated from in code documentation. (multiline comments under function headers) 

Date of generation: {datetime.now().strftime("%Y.%m.%d %H:%M:%S (YYYY.MM.DD HH:MM:SS)")} 
'''

def print_help(fun, indent=0):
    """Returns the __doc__ of the function as a string. 
    The string is formatted to be a markdown code block, and indented by the given number of tabs."""

    indent = indent * '\t'
    if not hasattr(fun, '__doc__') or fun.__doc__ is None:
        return ''
    text = fun.__doc__
    lines = text.split('\n')
    if text.startswith('_'):
        return ''
    if len(lines) < 2:
        return indent + '~~~\n' + indent + text.lstrip() + '\n' + indent + '~~~\n\n'
    spaces = len(lines[1]) - len(lines[1].lstrip())
    if lines[1].lstrip().startswith('*'):
        spaces -= 1

    text = indent + '~~~\n'
    for i, line in enumerate(lines):
        if (i == 0 or i == len(lines)-1) and line.lstrip() == '':
            continue
        line = line.removeprefix(spaces*' ')
        text += indent + line + '\n'
    
    return text + indent + '~~~\n\n'

def print_variable(name, value, use_type=False):
    """Returns the name and value of a variable as a formatted string."""

    var_type = str(getattr(value, '__class__', '\'None\'')).split("'")[-2]
    if len(str(value).replace('\n', '\\n')) < 80:
        text = str(value)
    else:
        text = str(value)[:80].replace('\n', ', ') + '...'
    if name == 'self' or var_type == 'NoneType':
        return name
    if type(value) == str:
        text = "'" + text + "'"
    if use_type:
        return f'**{name}**: *{var_type}* = `{text}`'
    return  f'{name}={text}'

def print_function(fun, indent=0):
    """Returns the function declaration as a string. 
    The string is formatted to be a markdown code block, and indented by the given number of tabs."""

    indent = indent * '\t'
    function_name = getattr(fun, '__name__', 'None')
    text = f'{function_name}('
    if type(fun) == type:
        fun = fun.__init__
        text = 'class ' + text
    else:
        text = 'def ' + text
    defaults = getattr(fun, '__defaults__', None)

    if defaults is not None:
        defaults = list(defaults)
    else:
        defaults = []

    if hasattr(fun, '__code__'):
            args = fun.__code__.co_varnames[:fun.__code__.co_argcount]
    else:
        args = []
    defaults = [None]*(len(args)-len(defaults)) + defaults
    for arg, default in zip(args, defaults):
        text += print_variable(arg, default) + ', '
    if len(args) > 0:
        text = text[:-2]
    return indent + '~~~Python\n' + indent + text + ')\n' + indent + '~~~\n\n'

def list_all(module, show_private=False):
    """Returns a list of all modules, classes, functions and variables in the given module. 
    Can hide private elements. (starting with '_') 
    Submodules are marked with a '#' in front of their name. """

    modules = []
    classes = []
    functions = []
    variables = []
    for i, element_name in enumerate(dir(module)):
        if element_name.startswith('__'):   # attribute
            continue
        if not show_private and element_name.startswith('_'):   # private
            continue
        element = getattr(module, element_name)

        if hasattr(element, '__package__'): # module (os, numpy, etc.)
            if hasattr(element, '__file__') and os.path.dirname(element.__file__) == os.path.dirname(module.__file__):
                modules.append('#' + element.__name__)
            else:
                modules.append(element.__name__)
            continue
        if hasattr(element, '__module__'):  # class or function
            if module.__name__ != element.__module__:   # different module
                if not element.__module__.split('.')[0] == module.__package__:    # not the same package
                    continue

            if element.__class__ == type:   # class
                classes.append(element_name)
                continue
            else:   # function
                functions.append(element_name)
                continue

        else:   # variable or other
            if not callable(element):
                variables.append(element_name)
                continue
            else:   # other
                functions.append(element_name)
                continue

    return modules, classes, functions, variables

def main():
    """Generates the documentation for the project."""

# Table of content
    cwd = os.getcwd()
    if 'Bubble_dynamics_simulation' not in cwd:
        file = open('Bubble_dynamics_simulation/DOCUMENTATION.md', 'w')
    else:
        file = open('DOCUMENTATION.md', 'w')
    file.write(header)
    file.write('## Table of contents\n\n')
    for i, module in enumerate(to_document):
        file.write(f'{i}. [**{module.__name__}.py**](#bookmark_{module.__name__})\n')

# Documenting
    for i, module in enumerate(to_document):
        file.write(f'<a name="bookmark_{module.__name__}"/>\n\n')
        file.write(f'## {module.__name__}.py\n\n')
        file.write(f'[Link to file]({module.__name__}.py)\n\n')
        file.write(f'{print_help(module, indent=0)}')

        modules, classes, functions, variables = list_all(module, document_private)

    # MODULES
        internal_modules = [mod for mod in modules if '#' in mod]
        external_modules = [mod for mod in modules if '#' not in mod]
        if len(internal_modules) > 0:
            file.write('Internal dependencies: \n')
        for mod in internal_modules:
            file.write(f' * {mod.replace("#", "")}.py\n')
        if len(external_modules) > 0:
            file.write('\nExternal dependencies: \n')
        for mod in external_modules:
            file.write(f' * {mod}\n')

    # CLASSES
        if len(classes) > 0:
            file.write('\n### Classes\n\n')
        for i, class_name in enumerate(classes):
            element = getattr(module, class_name)
            file.write(f'* **{element.__name__}**\n\n')
            file.write(print_function(element, indent=1))
            if hasattr(element, '__doc__') and element.__doc__ is not None:
                file.write(print_help(element, indent=1))
            for sub_element_name in dir(element):
                if class_name == 'dotdict':
                    continue
                if sub_element_name.startswith('__'):
                    continue
                if not document_private and sub_element_name.startswith('_'):
                    continue
                sub_element = getattr(element, sub_element_name)
                file.write(f'\t* *{sub_element.__name__}*\n\n')
                file.write(print_function(sub_element, indent=2))
                if help and hasattr(sub_element, '__doc__') and sub_element.__doc__ is not None:
                    file.write(print_help(sub_element, indent=2))   
    
    # FUNCTIONS
        if len(functions) > 0:
            file.write('\n### Functions\n\n')
        for i, function_name in enumerate(functions):
            function = getattr(module, function_name)
            file.write(f'* **{function.__name__}**\n\n')
            file.write(print_function(function, indent=1))
            if hasattr(function, '__doc__') and function.__doc__ is not None:
                file.write(print_help(function, indent=1))

    # VARIABLES
        if len(variables) > 0:
            file.write('\n### Global variables\n\nActual values might differ\n')
        for i, variable_name in enumerate(variables):
            variable = getattr(module, variable_name)
            file.write(' * ' + print_variable(variable_name, variable, use_type=True) + '\n')

    file.close()

if __name__ == '__main__':
    main()