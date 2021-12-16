#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
Created on    : Fri Dec 18 09:59:15 2020
Author        : Alejandro Martínez León
Mail          : [alejandro.marrtinezleon@uni-saarland.de, ale94mleon@gmail.com]
Affiliation   : Jochen Hub's Biophysics Group
Affiliation   : Faculty of NS, University of Saarland, Saarbrücken, Germany
===============================================================================
DESCRIPTION   :
DEPENDENCIES  :
===============================================================================
"""


from glob import glob
from rmsd import kabsch_rmsd
import pandas as pd
import numpy as np

#system('touch temp1.txt temp2.txt')
#system('sleep 5')
#system('rm temp1.txt')

def out_read(out):
    '''
    

    Parameters
    ----------
    out : str
        psi4 put file name.

    Returns
    -------
    check_freq, Gibbs free energy, exyz, xyz2RMSD_H.

    '''

    
    with open(out, 'rt', encoding='latin-1') as file:
        lines = file.readlines()

    check_freq = True
    G = 0    
    exyz = []
    xyz2RMSD_H =[]

    for i in range(len(lines)):
        if 'Final (previous) structure:' in lines[i]:
            for j in range(i+2,len(lines)):
                if 'Saving final (previous) structure.' in lines[j]: break
                split = lines[j].split()    
                exyz.append([split[0], float(split[1]), float(split[2]), float(split[3])])
                if 'H' not in lines[j]:
                    xyz2RMSD_H.append([float(split[1]), float(split[2]), float(split[3])])
        if 'Freq [cm^-1]' in lines[i]:
            freqs = lines[i].split()[2:]
            for f in freqs:
                '''
                 psi4 represent the imaginary number as 5i, 
                 if an error ocurreduring the float convertion, 
                 or it was converted but if less than cero (another error on sqrt) 
                 then frequ_chek =False
                '''
                try:
                    np.sqrt(float(f))
                except:
                    check_freq = False
        if 'Total G, Free enthalpy at  298.15 [K]' in lines[i]:
            G = float(lines[i].split('[K]')[1].split()[0])
        
            
            
            
    exyz = pd.DataFrame(exyz)
    xyz2RMSD_H = np.array(xyz2RMSD_H, dtype=float)    
    
    return check_freq, G, exyz, xyz2RMSD_H

#==============================================================================
#                                   Main program
#==============================================================================
outs = glob('*out')

if len(outs) == 0:
    exit('No se pueden encontrar los archivos .out')

first_names = {}
for out in outs:
    first_names.add(out.split('_Cell')[0])
first_names = list(first_names)

wrong_freq = []
for first_name in first_names:
    to_work = []
    for out in outs:
        if out.split('_Cell')[0] == first_name:
            to_work.append(out)  

    freq_checkes = []
    Gibbs_free_energies = []
    coords = []
    coord2RMSD_Hs =[]
    good = []

    for item in to_work:
        freq_check, Gibbs_free_energy, coord, coord2RMSD_H = (out_read(item))
        if freq_check==True:
            good.append(item)
            Gibbs_free_energies.append(Gibbs_free_energy)
            coords.append(coord)
            coord2RMSD_Hs.append(coord2RMSD_H)
        else:
            wrong_freq.append(item)
    paired = list(zip(good, Gibbs_free_energies, coords, coord2RMSD_Hs)) 
    ORDERED = sorted(paired, key=lambda x: x[1])
    to_print = []
    for i in range(len(ORDERED)):
        to_print.append([ORDERED[i][0], ORDERED[i][1]])
    to_print = pd.DataFrame(to_print, columns = ['opt_Molecule_from file', 'Gibbs (hartree/partícula)'])    
    print(to_print)
    print('\n')
    with open(ORDERED[0][0].split('_Cell')[0]+'.orden', 'wt') as final:
        to_print.to_string(final)

   
    rmsd_matrix = np.zeros((len(ORDERED),len(ORDERED)))
    for i in range(len(rmsd_matrix)):
        for j in range(len(rmsd_matrix)):
            if j < i:
                rmsd_matrix[i][j] = rmsd_matrix[j][i]
            elif j == i:
                rmsd_matrix[i][j] == 0.0
            else:
                # =============================================================
                #     Calculando RMSD
                # =============================================================
                P = ORDERED[i][3]
                Q = ORDERED[j][3]
                rmsd_matrix[i][j] = kabsch_rmsd(P, Q, translate = True)
    to_print_rmsd_matrix = pd.DataFrame(rmsd_matrix, index = list(good), columns = list(good))
    with open(ORDERED[0][0].split('_Cell')[0]+'.rmsd', 'wt') as final:
        to_print_rmsd_matrix.to_string(final) 
#========================Se exporta el .xyz==================
    with open(ORDERED[0][0].split('.')[0]+'.xyz', 'wt') as final:
        final.write(str(len(ORDERED[0][2]))+'\n\n')
        ORDERED[0][2].to_string(final, header=False, index=False)
 
#======================================================================
if len(wrong_freq) != 0:
    print('Las cálculos: %s presentaron frecuencias negativas o imaginarias.' % wrong_freq)
else:    
    print('Todas las frecuencias de los archivos en la carpeta de trabajo son mayores que cero.')                
                    
                
