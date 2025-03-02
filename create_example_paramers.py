"""This script cleans all python file in the "INP file examples" directory
and runst inp_data_extractor.extract() on all .inp files."""

import os
import importlib   # for reloading your own files
from termcolor import colored   # for colored error messages

# import inp_data_extractor.py as inp:
try:
    import inp_data_extractor as inp
except ImportError:
    try:
        from Bubble_dynamics_simulation import inp_data_extractor as inp
    except ImportError as error:
        print(colored(f'Error, \'inp_data_extractor.py\' not found', 'red'))
        raise error
    except Exception as error:
        print(colored(f'Error, \'inp_data_extractor.py\' failed to load', 'red'))
        raise error
except Exception as error:
    print(colored(f'Error, \'inp_data_extractor.py\' failed to load', 'red'))
    raise error
importlib.reload(inp)   # reload changes you made


# get directories
example_dir = os.path.join(os.getcwd(), 'INP file examples')
working_dir = os.getcwd()

# clean all python files in the directory
for file in os.listdir(example_dir):
    if file.endswith('.py'):
        if not file.endswith('create_example_parameers.py'):
            file_path = os.path.join(example_dir, file)
            os.remove(file_path)
            print(f'{file} removed')

# run inp_data_extractor.py on all .inp files in the directory
for file in os.listdir(example_dir):
    if file.endswith('.inp'):
        file_path = os.path.join(example_dir, file)
        base_name = os.path.splitext(file)[0]
        new_file_name = f'parameters_{base_name}.py'

        inp.extract(file_path)
        os.rename(f'parameters.py', new_file_name)
        os.replace(new_file_name, os.path.join(example_dir, new_file_name))
        print(f'{new_file_name} created')

inp.extract(os.path.join(example_dir, 'chem_Otomo2018_without_O_FIXED_by_Cantera.inp'))