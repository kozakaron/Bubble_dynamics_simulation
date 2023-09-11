"""________________________________Libraries________________________________"""

import numpy as np
from termcolor import colored
try:
    import data
    
except:
    print(print(colored(f'Error, \'data.py\' not found', 'red')))
comment = '!'

"""________________________________Functions________________________________"""

# Splits string by separator, and removes empty elements from the resulting list
def separate(string: str, separator: str=' ') -> list:
    ret = string.split(separator)
    ret.append('')
    while len(ret)>0 and ret[-1] == '':
        ret.remove('')
    return ret

# Finds the first line containing looking_for
def find(looking_for: str, lines: list) -> int:
    i = 0
    while i < len(lines) and not looking_for in lines[i]:
        i += 1
    if i >= len(lines):
        print(colored(f'Error in find(), \'{looking_for}\' not found', 'red'))
        return -1
    return i

# Changes order of lines in array, so that it will be ordered as new instead of original
def rearrange(array: list, original: list, new: list) -> list:
    newArray = []
    if len(array) != len(original) or len(array) != len(new):
        print(colored(f'Warning in rearrange(), lists have different lengths', 'yellow'))
    for orig in original:
        if not orig in new:
            print(colored(f'Warning in rearrange(), \'{orig}\' is in original, but not in new', 'yellow'))
    for i in range(len(new)):
        index = 0
        if new[i] in original:
            index = original.index(new[i])
        else:
            print(colored(f'Warning in rearrange(), \'{new[i]}\' is in new, but not in original', 'yellow'))
        newArray.append(array[index])
    return newArray

# Prints a 1D or 2D array or list
# Can leave comments above columns and after lines
def print_array(array, width=0, comments=[], columns=[], max_len=0):
    # empty array
    if len(array) == 0:
        return '[]'
    
    if type(array) == np.ndarray:
        array = list(array)
    if type(array[0]) == np.ndarray:
        for i, x in enumerate(array):
            array[i] = list(x)

    separator = ','
    arr_opener = '['
    arr_closer = ']'
    line_opener = '['
    line_closer = '],'
    remark = '#'

    text = ''
    
    if isinstance(array[0], list):
    # 2D array
        text += arr_opener + '\n'
        if columns != []:
            text += f'\t{remark}'
            for col in columns:
                text += f'{col: >{width}} '
            text += '\n'
        for i, x in enumerate(array):
            text += '\t' + line_opener
            for y in x:
                if type(y) == str:
                    y = '\'' + y + '\''
                text += f'{y: >{width}}' + separator
            if i != len(array)-1:
                text = text[:-1]
                text += line_closer
            else:
                if array[0] != []: text = text[:-1]
                text += line_closer
                text = text[:-1]
                text += ' '
            if len(comments) == len(array):
                text += f'    {remark} {comments[i]}'
            text += '\n'
        text += arr_closer
    # 1D array in 1 line
    elif max_len==0:
        if len(comments) == len(array):
            text += remark
            for i, comment in enumerate(comments):
                text += f'{comment: >{width}} '
            text += '\n'
        text += arr_opener
        for i, x in enumerate(array):
            if type(x) == str:
                x = '\'' + x + '\''
            text += f'{x: >{width}}' + separator
        text = text[:-1]
        text += arr_closer
    # 1D array in multiple lines
    else:
        text += arr_opener + '\n\t'
        i = 0
        while i < len(array):
            if len(comments) == len(array):
                text += remark
                for j in range(i, i+max_len):
                    if j >= len(array): break
                    text += f'{comments[j]: >{width}} '
                text += '\n\t '
            for j in range(i, i+max_len):
                if j >= len(array): break
                x = array[j]
                if type(x) == str:
                    x = '\'' + x + '\''
                text += f'{x: >{width}}' + separator
            
            i += max_len
            if i < len(array):
                text += '\n\t'
            else:
                text = text[:-1]
            
        text += '\n' + arr_closer
        
    return text
        
"""________________________________Line seperation, remove comments________________________________"""

def get_lines(text):  
    text = text.upper()
    text = text.replace('\\\\', '\\')
    text = text.replace('\t', ' ')
    text = text.replace('\r', ' ')
    text = text.replace('\v', '')
    text = text.replace('\b', '')
    text = text.replace('\f', '')
    text = text.replace('\a', '')
    text = text.replace('END', '\nEND')

    all_lines = separate(text, '\n')
    if len(all_lines) < 1:
        print(colored(f'Error, can not separate lines. Check if line ends are \'\\n\'', 'red'))
    all_lines.append('END')
    lines = []
    for line in all_lines:
        if line[0]==comment:
            continue
        if comment in line:
            line = line[:line.find(comment)]
        if line.replace(' ', '') != '':
            lines.append(line)
            
    return lines

"""________________________________Elements________________________________"""

def get_elements(lines):
    i = find('ELEM', lines)
    elements = []
    while not 'END' in lines[i]:
        line = lines[i]
        line = line.replace('ELEMENTS', '')
        line = line.replace('ELEM', '')
        line = line.replace('/', ' ')
        all_elements = separate(line, ' ')
        for element in all_elements:
            if element in data.W:
                if element in elements:
                    print(colored(f'Warning, element \'{element}\' is duplacated', 'yellow'))
                else:
                    elements.append(element)
            else:
                if element == '/': continue
                if element.replace('/', '').replace('+', '', 2).replace('-', '', 2).replace('.', '', 1).replace('E', '', 1).isnumeric(): continue
                print(colored(f'Warning, element \'{element}\' is not recognised, it won\'t be, included', 'yellow')) 
        i += 1
    
    return elements

"""________________________________Species data________________________________"""

def get_species(lines, elements):
  # Get species
    i = find('SPEC', lines)
    species = []
    while not 'END' in lines[i]:
        line = lines[i]
        line = line.replace('SPECIES', '')
        line = line.replace('SPEC', '')
        species += separate(line, ' ')
        i += 1

# Get W and lambda for species
    W = []
    lambdas = []
    components = np.zeros((len(species), len(elements)), dtype=np.int32)
    digits = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
    for i, specie in enumerate(species):
        if species.count(specie) > 1:
            print(colored(f'Error, specie \'{specie}\' is duplicated', 'red'))
        if specie in data.lambdas:
            lambdas.append(round(data.lambdas[specie], 5))
        else:
            lambdas.append(0.0)
            print(colored(f'Note, lambda value for specie \'{specie}\' is not in data.py: 0.0 is used', 'blue'))

        specie += '  '
        specie = specie.replace('EX', '')
        while len(specie) > 0 and specie[0] != ' ':
            if specie[0:2] in elements:
                if specie[2] in digits:
                    components[i][elements.index(specie[0:2])] += int(specie[2])
                    specie = specie[3:]
                else:
                    components[i][elements.index(specie[0:2])] += 1
                    specie = specie[2:]
            elif specie[0] in elements:
                if specie[1] in digits:
                    components[i][elements.index(specie[0])] += int(specie[1])
                    specie = specie[2:]
                else:
                    components[i][elements.index(specie[0])] += 1
                    specie = specie[1:]
            else:
                print(colored(f'Warning, {specie[0]} is not recognised in specie \'{species[i]}\', it will be neglected', 'yellow'))
                specie = specie[2:]

        w = 0.0
        for element, num in zip(elements, components[i]):
            w += num * data.W[element]
        W.append(round(w, 5))
        
    return species, W, lambdas

"""________________________________Thermodynamic data________________________________"""

def get_thermo(lines, species):
    TempRange = []
    a_low = []
    a_high = []
    materials = []

    i = find('THER', lines)
    if 'THERMO ALL' in lines[i]:
        i += 2

    while not 'END' in lines[i]:
        first_line = separate(lines[i], ' ')
        if 'TEMP' in first_line:
            i += 1
            while separate(lines[i], ' ')[0].replace('+', '', 2).replace('-', '', 2).replace('.', '', 1).replace('E', '', 1).isnumeric():
                i+=1        

        if not first_line[0] in materials:
            materials.append(first_line[0])
            if first_line[-2] == '0':
                first_line = first_line[:-2]
            else:
                first_line = first_line[:-1]
            TempRange.append( [float(first_line[-3]), float(first_line[-2]), float(first_line[-1])] )

            other_lines = [] 
            for j in range(1, 4): 
                line = lines[i+j]
                line = line[:-2]
                line = line.replace('E-', '_')
                line = line.replace('-', ' -')
                line = line.replace('_', 'E-')
                other_lines += separate(line, ' ')

            array_line = [float(num) for num in other_lines]
            a_high.append(array_line[0:7])
            a_low.append(array_line[7:14])

        i += 4

    TempRange = rearrange(TempRange, materials, species)
    a_low = rearrange(a_low, materials, species)
    a_high = rearrange(a_high, materials, species)

    return TempRange, a_low, a_high

"""________________________________Reactions________________________________"""

def get_reactions(lines, species):
  # Declare lists
    reactions = []
    numbers = []
    A = []
    B = []
    E = []

    ThirdBodyIndexes = []
    alfa = []

    PressureDependentIndexes = []
    LindemannIndexes = []
    ReacConst = []
    TroeIndexes = []
    Troe = []
    SRIIndexes = []
    SRI = []

    PlogIndexes = []
    Plog = []

  # Get reaction datas
    keywords = ['LOW', 'TROE', 'SRI', 'HIGH', 'REV', 'DUP', 'LT', 'TDEP', 'XSMI', 'PLOG', 'FORD', 'RORD', 'MOME', 'EXCI', 'JAN', '/']
    i = find('REAC', lines) + 1

    while not 'END' in lines[i]: 
        reaction_line = lines[i]
        reaction_line = reaction_line.replace('<=>', '=')
        reaction_line = reaction_line.replace('=>', '>')
        reaction_line = separate(reaction_line, ' ')
        reactions.append(''.join(reaction_line[:-3]).replace('>', '=>'))
        numbers.append(len(reactions)-1)
        A.append(float(reaction_line[-3]))
        B.append(float(reaction_line[-2]))
        E.append(float(reaction_line[-1]))
        
        line = lines[i]
        isThirdBody = isPressureDependent = isTroe = isSRI = isPLOG = False
        isManualThirdBodyCoefficients = False
        
        if '(+M)=' in line.replace(' ',''):
            isPressureDependent = True
            isThirdBody = True
        elif '+M=' in line.replace(' ',''): 
            isThirdBody = True
        
        if not any([keyword in lines[i] for keyword in keywords+['END']]):
            i += 1
            line = lines[i]
        
        while any([keyword in line for keyword in keywords+['END']]): #This cycle steps one reaction.
            if 'END' in line:
                break
            elif 'LOW' in line:
                isThirdBody = True
                isPressureDependent = True
                line = line.replace('LOW', '')
                line = line.replace('/', '')
                line = separate(line, ' ')
                ReacConst.append([float(line[-3]), float(line[-2]), float(line[-1])])
            elif 'TROE' in line:
                isTroe = True
                line = line.replace('TROE', '')
                line = line.replace('/', '')
                line = separate(line, ' ')
                if len(line) >= 4:
                    Troe.append([float(line[-4]), float(line[-3]), float(line[-2]), float(line[-1])])
                else:
                    Troe.append([float(line[-3]), float(line[-2]), float(line[-1]), 1e300]) # last coeff is inf
                    print(colored(f'Note, no T*** in line {i} (\'{lines[i]}\') in reaction \'{reactions[-1]}\': T***=inf is used', 'blue'))
            elif 'SRI' in line:
                isSRI = True
                line = line.replace('SRI', '')
                line = line.replace('/', '')
                line = separate(line, ' ')
                if len(line) == 5:
                    SRI.append([float(line[-5]), float(line[-4]), float(line[-3]), float(line[-2]), float(line[-1])])
                else:
                    SRI.append([float(line[-3]), float(line[-2]), float(line[-1]), 1.0, 0.0]) # d=1, e=0
                    print(colored(f'Note, no d or e in line {i} (\'{lines[i]}\') in reaction \'{reactions[-1]}\': d=1, e=0 is used', 'blue'))
            elif 'PLOG' in line:
                if not isPLOG:
                    PlogIndexes.append(numbers[-1])
                    isPLOG = True
                line = line.replace('PLOG', '')
                line = line.replace('MX', '')
                line = line.replace('SP', '')
                line = line.replace('/', '')
                line = separate(line, ' ')
            
                Plog.append([float(line[-4]), float(line[-3]), float(line[-2]), float(line[-1])])
            elif '/' in line:
                line = line.replace('/ ', ' ')
                line = line.replace('/', ' ')
                line = separate(line, ' ')
                alfa_line = np.ones((len(species)), dtype=np.float64)
                isManualThirdBodyCoefficients = True
                j = 0
                while j < len(line):
                    if line[j] in species:
                        alfa_line[species.index(line[j])] = round(float(line[j+1]), 5)
                    else:
                        print(colored(f'Warning, third body \'{line[j]}\' is not in species in line {i} (\'{lines[i]}\') in reaction \'{reactions[-1]}\'', 'yellow'))
                    j += 2
                alfa.append(alfa_line)
            elif 'DUP' in line:
                line = line.replace('DUP','')
            else:
                keyword = keywords[[keyword in lines[i] for keyword in keywords].index(True)]
                print(colored(f'Warning, keyword \'{keyword}\' is not supported in line {i} (\'{line}\')', 'yellow'))
            if not any([keyword in line for keyword in keywords+['END']]):
                i += 1
                line = lines[i]

        if isPressureDependent:
            PressureDependentIndexes.append(numbers[-1])
            if isTroe:
                TroeIndexes.append(numbers[-1])
            elif isSRI:
                SRIIndexes.append(numbers[-1])
            else:
                LindemannIndexes.append(numbers[-1])
        if isThirdBody:
            ThirdBodyIndexes.append(numbers[-1])
            if not isManualThirdBodyCoefficients:
                alfa_line = np.ones((len(species)), dtype=np.float64)
                alfa.append(alfa_line)
        if isPLOG and len(Plog) % 3 != 0:
            print(colored(f'Warning, only 3 lines of PLOG is supported in reaction \'{reactions[-1]}\'', 'yellow'))

        
    try:
        Plog = np.array(Plog)
    except:
        Plog = [[]]

    if ReacConst == []: ReacConst = [[]]
    if Troe == []: Troe = [[]]
    if SRI == []: SRI = [[]]

    A = np.array(A)
    B = np.array(B)
    E = np.array(E)
    ReacConst = np.array(ReacConst)
    
  # Correct units
    i = find('REAC', lines)
    if 'MOLEC' in lines[i]: # MOLECULE
        print(colored(f'Note, pre-exponential factor is modified from units of [cm^3/molecule/s] to [cm^3/mol/s]', 'blue'))
        A /= 6.02214e23
        if len(PressureDependentIndexes) != 0: ReacConst[:, 0] /= 6.02214e23
        Plog[:, 1]  /= 6.02214e23
        # Avogadro's number: N_A = 6.02214e23 [-]
    elif 'MOL' in lines[i]:
        print(colored(f'Note, pre-exponential factor is modified from units of [cm^3/mol/s] to [cm^3/mol/s]', 'blue'))
    if 'CAL' in lines[i]:
        print(colored(f'Note, activation energy (E_i) is modified from units of [cal/mol] to [cal/mol]', 'blue'))
    elif 'KCAL' in lines[i]:
        print(colored(f'Note, activation energy (E_i) is modified from units of [kcal/mol] to [cal/mol]', 'blue'))
        E /= 1000.0 # [kcal/mol -> cal/mol]
        if len(PressureDependentIndexes) != 0: ReacConst[:, 2] /= 1000.0
        if len(PlogIndexes) != 0: 
            Plog[:, 3]  /= 1000.0
    elif 'JOU' in lines[i]:
        print(colored(f'Note, activation energy (E_i) is modified from units of [J/mol] to [cal/mol]', 'blue'))
        E /= 4.184 # [J/mol -> cal/mol]
        if len(PressureDependentIndexes) != 0: ReacConst[:, 2] /= 4.184
        if len(PlogIndexes) != 0: 
            Plog[:, 3]  /= 4.184
    elif 'KJOU' in lines[i]:
        print(colored(f'Note, activation energy (E_i) is modified from units of [kJ/mol] to [cal/mol]', 'blue'))
        E *= 1000.0 # [kJ/mol -> J/mol]
        E /= 4.184 # [J/mol -> cal/mol]
        if len(PressureDependentIndexes) != 0: ReacConst[:, 2] *= 1000.0 / 4.184
        if len(PlogIndexes) != 0: 
            Plog[:, 3]  *= 1000.0 / 4.184
    elif 'KELV' in lines[i]:
        print(colored(f'Note, activation energy (E_i) is modified from units of [K] to [cal/mol]', 'blue'))
        E *= 1.9872 # [K -> cal/mol]
        if len(PressureDependentIndexes) != 0: ReacConst[:, 2] *= 1.9872
        if len(PlogIndexes) != 0: 
            Plog[:, 3]  *= 1.9872
        # Universal gas constant: impal = 1.987 [cal/mol/K]
    elif 'EVOL' in lines[i]:
        print(colored(f'Note, activation energy (E_i) is modified from units of [eV/mol] to [cal/mol]', 'blue'))
        E *= 1.602176634e-19 # [eV/mol -> J/mol]
        E /= 4.184 # [J/mol -> cal/mol]
        if len(PressureDependentIndexes) != 0: ReacConst[:, 2] *= 1.602176634e-19 / 4.184
        if len(PlogIndexes) != 0: 
            Plog[:, 3]  *= 1.602176634e-19 / 4.184

    for i in range(len(E)):
        E[i] = round(E[i], 5)
    
    return (reactions, A, B, E,
            ThirdBodyIndexes, alfa,
            PressureDependentIndexes,
                LindemannIndexes, ReacConst,
                TroeIndexes, Troe,
                SRIIndexes, SRI,
            PlogIndexes, Plog)

"""________________________________Reaction matrixes (nu)________________________________"""

def get_nu(reactions, species, W):
    nu_forward = np.zeros((len(reactions), len(species)), dtype=int)
    nu_backward = np.zeros((len(reactions), len(species)), dtype=int)
    IrreversibleIndexes = []
    digits = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']

    for x, reaction in enumerate(reactions):
        if '>' in reaction:
            IrreversibleIndexes.append(x)
        reaction = reaction.replace('+M', '')
        reaction = reaction.replace('(+M)', '')
        reaction = reaction.replace('()', '')
        reaction = reaction.replace('>', '=')
        forward = separate(reaction, '=')[0]
        backward = separate(reaction, '=')[1]
        forward = separate(forward, '+')
        backward = separate(backward, '+')
        for f in forward:
            num = 1
            if f[0] in digits:
                num = int(f[0])
                f = f[1:]
            if f in species:
                nu_forward[x][species.index(f)] += num
            else:
                if f =='E':
                    print(colored(f'Warning, electron (E) in reaction {x} (\'{reactions[x]}\') is ignored', 'yellow'))
                elif f == 'HV':
                    print(colored(f'Warning, photon (HV) in reaction {x} (\'{reactions[x]}\') is ignored', 'yellow'))
                else:
                    print(colored(f'Warning, \'{f}\' in reaction {x} (\'{reactions[x]}\') is not in species, the reaction is ignored', 'yellow'))

        for b in backward:
            num = 1
            if b[0] in digits:
                num = int(b[0])
                b = b[1:]
            if b in species:
                nu_backward[x][species.index(b)] += num
            else:
                if f =='E':
                    print(colored(f'Warning, electron (E) in reaction {x} (\'{reactions[x]}\') is ignored', 'yellow'))
                elif f == 'HV':
                    print(colored(f'Warning, photon (HV) in reaction {x} (\'{reactions[x]}\') is ignored', 'yellow'))
                else:
                    print(colored(f'Warning, \'{f}\' in reaction {x} (\'{reactions[x]}\') is not in species, the reaction is ignored', 'yellow'))

    nu = nu_backward - nu_forward
    for i in range(0, len(reactions)):
        if abs(sum(nu[i] * W)) > 1e-5:
            print(colored(f'Warning, nonconsistent reaction {i} (\'{reactions[i]}\')', 'yellow'))
    
    return nu_forward, nu_backward, IrreversibleIndexes

"""________________________________Printing to file________________________________"""

def extract(path):
  # Open file
    try:
        file = open(path, 'r')
        text = file.read()
        model = path.split('\\')[-1][:-4]
    except:
        print(colored(f'Error, \'{path}\' not found', 'red'))
    
  # Extract data
    lines = get_lines(text)
    elements = get_elements(lines)
    species, W, lambdas = get_species(lines, elements)
    TempRange, a_low, a_high = get_thermo(lines, species)
    (reactions, A, B, E,
        ThirdBodyIndexes, alfa,
        PressureDependentIndexes,
            LindemannIndexes, ReacConst,
            TroeIndexes, Troe,
            SRIIndexes, SRI,
        PlogIndexes, Plog) = get_reactions(lines, species)
    nu_forward, nu_backward, IrreversibleIndexes = get_nu(reactions, species, W)
    
  # Create parameters.py
    line_start = '\n\"\"\"________________________________'
    line_end = '________________________________\"\"\"\n'
    text = ''
    
    # Physical constants, Elements and Species data
    text += data.header + '\n\n'
    text += f'model = \'{model}\'\n'
    text += f'import numpy as np\n'
    text += line_start + 'Physical constants' + line_end
    text += data.physical_constants + '\n\n'
    
    # Species and elements
    text += line_start + 'Species' + line_end
    text += f'elements = np.array({print_array(elements, 0)})\n'
    if len(species) < 16:
        text += f'#                   {print_array([i for i in range(len(species))], 10)[1:-1]}\n'
        max_len = 0
    else:
        max_len = 10
    text += f'species = np.array({print_array(species, 10, max_len=max_len)})\n'
    text += f'# molar mass [g/mol]\n'
    text += f'W = np.array([      {print_array(W, 10, max_len=max_len)[1:]}, dtype=np.float64)\n'
    text += f'# thermal conductivity [W / m / K]\n'
    text += f'lambdas = np.array({print_array(lambdas, 10, max_len=max_len)}, dtype=np.float64)\n'
    text += f'index = dict(\n\t'
    for i, specie in enumerate(species):
        text += f'{specie: >6}={i: >2}'
        if (i+1) % 10 == 0 and i != len(species)-1:
            text += ',\n\t'
        elif i != len(species)-1:
            text += ', '
    text += '\n)\n'
    text += f'indexOfWater = ' + str(species.index('H2O')) + '\n'
    text += f'K = {len(W)}   # Number of species\n\n'
    
    # NASA polynomials
    text += line_start + 'NASA polynomials' + line_end
    text += f'N = 5    # degree of polynomials\n'
    text += f'TempRange = np.array('+ print_array(TempRange, 8, species, ['T_low', 'T_high', 'T_mid']) + ', dtype=np.float64)\n\n'
    text += f'# LOW NASA coefficients\n'
    text += f'a_low = np.array('+ print_array(a_low, 16, species, ['a_1', 'a_2', 'a_3', 'a_4', 'a_5', 'a_6', 'a_7']) + ', dtype=np.float64)\n\n'
    text += f'# LOW NASA coefficients\n'
    text += f'a_high = np.array('+ print_array(a_high, 16, species, ['a_1', 'a_2', 'a_3', 'a_4', 'a_5', 'a_6', 'a_7']) + ', dtype=np.float64)\n\n'
    
    # Reaction constants
    text += line_start + 'Reaction constants' + line_end
    text += f'I = {len(reactions)}    # Number of reactions\n'
    text += f'# Pre-exponential factors [cm^3/mol/s v 1/s]\n'
    text += f'A = np.array('+ print_array(A, 20, max_len=5) + ', dtype=np.float64)\n\n'
    text += f'# Temperature exponent [-]\n'
    text += f'b = np.array('+ print_array(B, 20, max_len=5) + ', dtype=np.float64)\n\n'
    text += f'# Activation energy [cal/mol]\n'
    text += f'E = np.array('+ print_array(E, 20, max_len=5) + ', dtype=np.float64)\n\n'
    
    # Reaction matrixes
    text += line_start + 'Reaction matrixes' + line_end
    text += f'# Forward reaction matrix\n'
    text += f'nu_forward = np.array('+ print_array(nu_forward, 4, [f'{x:>2}. {reaction}' for x, reaction in enumerate(reactions)], species) + ', dtype=np.float64)\n\n'
    text += f'# Backward reaction matrix\n'
    text += f'nu_backward = np.array('+ print_array(nu_backward, 4, [f'{x:>2}. {reaction}' for x, reaction in enumerate(reactions)], species) + ', dtype=np.float64)\n\n'
    text += f'nu = nu_backward - nu_forward\n\n'
    
    # Three-body reactions
    text += line_start + 'Three-body reactions' + line_end
    text += f'ThirdBodyIndexes = np.array({print_array(ThirdBodyIndexes, 4)}, dtype=np.int64)\n'
    text += f'ThirdBodyCount = {len(ThirdBodyIndexes)}\n\n'
    text += f'# third-body efficiency factors\n'
    text += f'alfa = np.array('+ print_array(alfa, 8, [f'{x:>2}. {reactions[x]}' for x in ThirdBodyIndexes], species) + ', dtype=np.float64)\n\n'
    
    # Irreversible reactions
    text += line_start + 'Irreversible reactions' + line_end
    text += f'IrreversibleIndexes = np.array({print_array(IrreversibleIndexes, 4)}, dtype=np.int64)\n'
    text += f'IrreversibleCount = {len(IrreversibleIndexes)}\n\n'
    
    # Pressure-dependent reactions
    text += line_start + 'Pressure-dependent reactions' + line_end
    text += f'PressureDependentIndexes = np.array({print_array(PressureDependentIndexes, 4)}, dtype=np.int64)\n'
    text += f'PressureDependentCount = {len(PressureDependentIndexes)}\n\n'
    text += f'LindemannIndexes = np.array({print_array(LindemannIndexes, 4)}, dtype=np.int64)\n'
    text += f'LindemannCount = {len(LindemannIndexes)}\n\n'
    text += f'# Fall-off parameters\n'
    text += f'ReacConst = np.array('+ print_array(ReacConst, 18, [f'{x:>2}. {reactions[x]}' for x in PressureDependentIndexes], ['A_0', 'b_0', 'E_0']) + ', dtype=np.float64)\n\n'
    
    text += f'TroeIndexes = np.array({print_array(TroeIndexes, 4)}, dtype=np.int64)\n'
    text += f'TroeCount = {len(TroeIndexes)}\n\n'
    text += f'# Troe parameters\n'
    text += f'Troe = np.array('+ print_array(Troe, 18, [f'{x:>2}. {reactions[x]}' for x in TroeIndexes], ['alfa',  'T***',  'T*',  'T**']) + ', dtype=np.float64)\n\n'
    
    text += f'SRIIndexes = np.array({print_array(SRIIndexes, 4)}, dtype=np.int64)\n'
    text += f'SRICount = {len(SRIIndexes)}\n\n'
    text += f'# SRI parameters\n'
    text += f'SRI = np.array('+ print_array(SRI, 18, [f'{x:>2}. {reactions[x]}' for x in SRIIndexes], ['a',  'b',  'c',  'd', 'e']) + ', dtype=np.float64)\n\n'
    
    text += f'PlogIndexes = np.array({print_array(PlogIndexes, 4)}, dtype=np.int64)\n'
    text += f'PlogCount = {len(PlogIndexes)}\n\n'
    text += f'# PLOG parameters\n'
    text2 = f'Plog = np.array([\n#'+print_array(['P_1',  'A_1',  'b_1',  'E_1'],12)[1:-1]+'\n'
    text2 = text2.replace("'", '')
    text2 = text2.replace(",", ' ')
    text += text2
    for i in range(len(PlogIndexes)-1):
        text += print_array(Plog[3*i], 10)+', '+str(f'#{PlogIndexes[i]:>2}. {reactions[PlogIndexes[i]]}')+'\n'
        text += print_array(Plog[3*i+1], 10)+', \n'
        text += print_array(Plog[3*i+2], 10)+', \n'
    if(len(PlogIndexes)>0):
        text += print_array(Plog[3*len(PlogIndexes)-3], 10)+','+str(f' #{PlogIndexes[len(PlogIndexes)-1]:>2}. {reactions[PlogIndexes[len(PlogIndexes)-1]]}')+ '\n'
        text += print_array(Plog[3*len(PlogIndexes)-2], 10)+', \n'
        text += print_array(Plog[3*len(PlogIndexes)-1], 10)
    text += '\n], dtype=np.float64)\n'
    
    text = text.replace('\t', '    ')
    file = open('parameters.py', 'w', encoding='utf8')
    file.write(text)
    file.close()
    
    print(f'model: {model}')
    print(f'File \'parameters.py\' succesfully created')