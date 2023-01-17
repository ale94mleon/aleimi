# -*- coding: utf-8 -*-
"""
===============================================================================
Created on    : Tue Dec  8 00:51:27 2020
Author        : Alejandro Martínez León
Mail          : [alejandro.marrtinezleon@uni-saarland.de, ale94mleon@gmail.com]
Affiliation   : Jochen Hub's Biophysics Group
Affiliation   : Faculty of NS, University of Saarland, Saarbrücken, Germany
===============================================================================
DESCRIPTION   :
DEPENDENCIES  :
===============================================================================
"""
# Se podria especificar los dos directorios en donde se desea comparar
# o un directorio de referencia (el de mayor nivel de teoria MP2)
# y luego todos los demas en los que se desee comparar con este, pues puedo tener
# B3LYP, HF y unos cuantos y entonces tendria que correr varias veces el mismo
# calcualo, seria recomendable que el directoria, especifique el nievel de teria 
# empleado para tomar esto como nombre.
# me queda implementar las funciones de los con signo y sin signo. Que lo que deberia es 
import calc_bond_angle as calc
import math as mt
import pandas as pd
from os import walk
import glob as glob

# =============================================================================
#      USER SPECIFICATIONS
# By defualt a single point energy is performed. You can specified if 
# freq is used or not.
# =============================================================================
higher_level_dir = ['MP2'] # reference
low_level_dir = ['HF']
check_function = calc.psi4_freq_check ## I can use psi4_opt_check insted of psi4_freq_check to compare the structures neve mind the imaginary modes

# When the list is used -1 need to be substrated. Because python array start at 0 not 1
bond_list = [(10,11),
             (11,12),
             (12,13),
             (13,14),
             (14,15),
             (15,10),
             (10,8),
             (8,9),
             (8,7),
             (7,5),
             (5,6),
             (5,2),
             (2,4),
             (2,3),
             (2,1),
             (1,16),
             (1,17),
             (1,18)
                     ]
angle_bond_list = [(11,10,8),
                   (10,8,9), 
                   (9,8,7),
                   (7,5,6),
                   (6,5,2),
                   (5,2,3),
                   (5,2,4),
                   (5,2,1),
                   (2,1,16),
                   (2,1,17),
                   (2,1,18)
                               ]
dihedral_list = [(11,10,8,7),
                  (10,8,7,5),
                  (8,7,5,2),
                  (7,5,6,19),
                  (7,5,2,1),
                  (5,2,1,16)
                             ]
# =============================================================================

py_bond_list = [list(map(lambda x:x-1, item)) for item in bond_list]
py_angle_bond_list = [list(map(lambda x:x-1, item)) for item in angle_bond_list]
py_dihedral_list = [list(map(lambda x:x-1, item)) for item in dihedral_list]



# Match each conformer/mols with the reference level
directories = higher_level_dir + low_level_dir
psi4_outs_alls = []
for directory in directories:
    psi4_outs_per_level = []
    for root, dirs, files in (walk('./'+directory)):    
        #print(root)
        #print(dirs)
        #print(files)
        outs = glob.glob(root+'/*.out')
        for out in outs:
            if 'myjob' not in out:
                psi4_outs_per_level.append(out)
    psi4_outs_alls.append(psi4_outs_per_level)

match = []
for item_higher in psi4_outs_alls[0]:
    match_temp = [item_higher]
    for lower_level in psi4_outs_alls[1:]:
        for item_lower in lower_level:
            if item_higher.split('/')[-1] == item_lower.split('/')[-1]:# check if the files have the same name
                match_temp.append(item_lower)
    if len(match_temp) > 1: match.append(match_temp) # if match something

# Each item in match contein a mol at different level of theories
# mol[0] has the higher level of theory, if there is some problem with the file
# delete the item of match
# is important to sorted in reverse for the pop function
idx_2_delete = sorted([i for i, mol in enumerate(match) if not check_function(mol[0])], reverse=True)

# Only save the level tested (the higher), the rest for these molecules will not tested
wrong_freq = [match.pop(idx)[0] for idx in idx_2_delete]



match_final = []
for mol in match:
    match_final_temp = [mol[0]]
    for lower_level in mol[1:]:
        if not check_function(lower_level):
            wrong_freq.append(lower_level)
        else:
            match_final_temp.append(lower_level)
    if len(match_final_temp) > 1: match_final.append(match_final_temp) # If more that only the higher level

with open('Wrong_freq.txt', 'wt') as file:
    file.write('The following files had freq or opt problems. For such molecules that the higher level of theory failed, it was not checked the rest level of theories.\n')
    [file.write(w+'\n') for w in wrong_freq]       



energy = []
RMSD_matrix = []
for mol in match_final:
    mol_energy = []
    mol_RMSD_matrix = []
    for level in mol:    
        mol_energy.append(calc.psi4_free_energy(level))
        mol_RMSD_matrix.append(calc.RMSD_xyz(level, mol[0]))
    energy.append(mol_energy)
    RMSD_matrix.append(mol_RMSD_matrix)

df_energy = pd.DataFrame(energy, index = [mol[0].split('/')[-1].split('.')[0] for mol in match_final], columns = [level.split('/')[1] for level in match_final[0]])
df_RMSD_matrix = pd.DataFrame(RMSD_matrix, index = [mol[0].split('/')[-1].split('.')[0] for mol in match_final], columns = [level.split('/')[1] for level in match_final[0]])
energy_RMSD_writer = pd.ExcelWriter('Energies_RMSD.ods')
df_energy.to_excel(energy_RMSD_writer, 'Energies')
df_RMSD_matrix.to_excel(energy_RMSD_writer, 'RMSD_xyz_ref_'+match_final[0][0].split('/')[1])
energy_RMSD_writer.save()


for i, mol in enumerate(match_final):   
    for j, level in enumerate(mol):
        xyz = calc.psi4out2coords(level)
        bonds = [calc.calc_bond_dist(*[xyz[item] for item in bond]) for bond in py_bond_list]
        angles = [calc.calc_angle_bond(*[xyz[item] for item in angle])*180/mt.pi for angle in py_angle_bond_list]
        dihedrals = [calc.calc_dihedral(*[xyz[item] for item in dihedral])*180/mt.pi for dihedral in py_dihedral_list]
        
        match_final[i][j] = [match_final[i][j], bonds, angles, dihedrals]
        


'''
Now, the structure of match_final is the following:
                                                            LevelS
    mol             MP2                                 HF                              B3LYP                       ...
    1    [name, bond, angle, dihedral]   [name, bond, angle, dihedral]          [name, bond, angle, dihedral]
    2    [name, bond, angle, dihedral]   [name, bond, angle, dihedral]          [name, bond, angle, dihedral]
    3    [name, bond, angle, dihedral]   [name, bond, angle, dihedral]          [name, bond, angle, dihedral]
    4   [name, bond, angle, dihedral]   [name, bond, angle, dihedral]          [name, bond, angle, dihedral]
    ...
    the next stepo want to reshape the data:
    for mol 1 
                                   Levels
    bonds    MP2                HF                                B3LYP
    1       x1                 x1                               x1
    2       x2                 x2                               x2
    3       x3                 x3                               x3
    4       x4                 x4                               x4
    ...
    and the same for the rest of the molecules. Each molecule will be exported to an excel file with three different sheets: bond, angles, dihedral
'''
for mol in match_final:
    columns = [item[0] for item in mol]
    # I need to rotate the data, this is done with transpouse ()
    df_bond = pd.DataFrame(calc.transpouse([item[1] for item in mol]), index = bond_list, columns = columns)
    df_angle = pd.DataFrame(calc.transpouse([item[2] for item in mol]), index = angle_bond_list, columns = columns)
    df_dihedral = pd.DataFrame(calc.transpouse([item[3] for item in mol]), index = dihedral_list, columns = columns)
    
    name = mol[0][0].split('/')[-1].split('.')[0]#+'benchmark'
    errors_bond = []
    errors_angle = []
    errors_dihedral = []
    writer = pd.ExcelWriter(name + '_benchmark.ods')
    for level in columns: #[1:] and them add a 0 in at the first position MP2/MP2

        dif_bond = list(df_bond[level] - df_bond[columns[0]]) # lower level - higher level(the higher level is constant for each molecule)
        dif_angle = list(df_angle[level] - df_angle[columns[0]])
        
        '''
        The following convention will be used.
        The dihedral angale is always in the interval -pi to pi.
        If both dihedrals are positive or negative at the same time:
            A negative difference means that the lower level move clockwise respect to the higher level
            And a positive difference means that the lower level move counterclockwise respect to the higher level.
        For mix signs:
        The correction down writed. The idea is take the same convention when both dihedrals are positive or negative at the same time
        This correction is needed because the real difference between 
        angels like -179 and 179 are not well estimated.
        The only differ in 2 degrees on the chemical structure.
        But -179-179 is not that.
        '''
        
        dif_dihedral = []
        for i in range(len(df_dihedral[level])): # Because I want to loop for each row
            lower_level = df_dihedral[level][i]
            higher_level = df_dihedral[columns[0]][i]
            if lower_level*higher_level < 0: # mix signs
                item1 = abs(higher_level) + abs(lower_level)
                item2 = 360 - item1
                if item1 < item2: # you need to get the small angle. Visualize the problem. Really important this check
                    if higher_level > 0:
                        item = -item1
                    else:
                        item = item1
                else:
                    if higher_level > 0:
                        item = item2
                    else:
                        item = -item2       
            else:
                item = lower_level - higher_level
            dif_dihedral.append(item)
        
        
        errors_bond.append(calc.errors(dif_bond))
        errors_angle.append(calc.errors(dif_angle))
        errors_dihedral.append(calc.errors(dif_dihedral))
    df_bond_errors = pd.DataFrame(calc.transpouse(errors_bond), index = ['RMSE', 'MSE', 'MUE'], columns=columns)
    df_angle_errors = pd.DataFrame(calc.transpouse(errors_angle), index = ['RMSE', 'MSE', 'MUE'], columns=columns)
    df_dihedral_errors = pd.DataFrame(calc.transpouse(errors_dihedral), index = ['RMSE', 'MSE', 'MUE'], columns=columns)
    
    
    df_bond_to_print = pd.concat([df_bond, df_bond_errors])
    df_angle_to_print = pd.concat([df_angle, df_angle_errors])
    df_dihedral_to_print = pd.concat([df_dihedral, df_dihedral_errors])

    df_bond_to_print.to_excel(writer, 'bond')
    df_angle_to_print.to_excel(writer, 'angle')
    df_dihedral_to_print.to_excel(writer, 'dihedral')
    writer.save()