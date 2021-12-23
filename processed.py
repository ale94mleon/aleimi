#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
Created on    : 2020-2023
Author        : Alejandro Martínez León
Mail          : [alejandro.martinezleon@uni-saarland.de, ale94mleon@gmail.com]
Affiliation   : Jochen Hub's Biophysics Group
Affiliation   : Faculty of NS, University of Saarland, Saarbrücken, Germany
===============================================================================
DESCRIPTION   :
DEPENDENCIES  :
===============================================================================
"""
import pandas as pd
import numpy as np
from glob import glob
import os
from rmsd import kabsch_rmsd
import tempfile
from aleimi import OBconvert, templates, tools

def psi4_out_read(out):
    """
    Parameters
    ----------
    out : str
        psi4 put file name.
    Returns
    -------
    check_freq, Gibbs free energy, exyz, xyz2RMSD_H.
    """
    
    with open(out, 'rt', encoding='latin-1') as file:
        lines = file.readlines()

    if lines[-1].strip() == '*** Psi4 exiting successfully. Buy a developer a beer!':
        check_end = True
    else:
        check_end  = False
    check_freq = False
    E = 0
    G = 0    
    exyz = []
    xyz2RMSD_H =[]

    for i in range(len(lines)):
        if 'Final (previous) structure:' in lines[i]:
            E = float(lines[i-1].split()[-1])
            for j in range(i+2,len(lines)):
                if 'Saving final (previous) structure.' in lines[j]: break
                split = lines[j].split()    
                exyz.append([split[0], float(split[1]), float(split[2]), float(split[3])])
                if 'H' not in lines[j]:
                    xyz2RMSD_H.append([float(split[1]), float(split[2]), float(split[3])])
        if 'Freq [cm^-1]' in lines[i]:
            check_freq = True
            freqs = lines[i].split()[2:]
            for f in freqs:
                """
                 psi4 represent the imaginary number as 5i, 
                 if an error ocurreduring the float convertion, 
                 or it was converted but if less than cero (another error on sqrt) 
                 then check_freq = False
                """
                try:
                    np.sqrt(float(f))
                except:
                    check_freq = False
                    break
        if 'Total G, Free enthalpy at  298.15 [K]' in lines[i]:
            G = float(lines[i].split('[K]')[1].split()[0])
        
    exyz = pd.DataFrame(exyz)
    xyz2RMSD_H = np.array(xyz2RMSD_H, dtype=float)   
    return check_end, check_freq, E, G, exyz, xyz2RMSD_H

def main(SubDirs = True, engine = 'psi4', xyz_out = False, parameterize_path = './parameterize'):

    if SubDirs:
        outs = [out for out in glob('*/*.out') if 'myjob' not in out]

        first_names = set()
        for out in outs:
            first_names.add(os.path.basename(out).split('_Cell')[0])
        first_names = list(first_names)
    else:
        outs = [out for out in glob('*.out') if 'myjob' not in out]
        
        first_names = set()
        for out in outs:
            first_names.add(out.split('_Cell')[0])
        first_names = list(first_names)
    if not outs: raise FileNotFoundError('No se pueden encontrar los archivos .out')

    if engine == 'psi4':
        reader = psi4_out_read
    elif engine == 'orca':
        pass
    elif engine == 'gaussian':
        pass
    else:
        pass

    wrong_freq = []
    wrong_end = []

    coord_lower_energy = {}
    for first_name in first_names:
        to_work = []
        if SubDirs:
            for out in outs:
                if os.path.basename(out).split('_Cell')[0] == first_name:
                    to_work.append(out)
        else:
            for out in outs:
                if out.split('_Cell')[0] == first_name:
                    to_work.append(out)  

        Gibbs_free_energies = []
        coords = []
        coord2RMSD_Hs =[]
        good = []

        for item in to_work:
            check_end, check_freq, E_total_energy, Gibbs_free_energy, coord, coord2RMSD_H = reader(item)
            if check_end:
                if check_freq==True:
                    good.append(os.path.basename(item))
                    Gibbs_free_energies.append(Gibbs_free_energy)
                    coords.append(coord)
                    coord2RMSD_Hs.append(coord2RMSD_H)
                else:
                    wrong_freq.append(item)
            else:
                wrong_end.append(item)
        if good:
            paired = list(zip(good, Gibbs_free_energies, coords, coord2RMSD_Hs)) 
            ORDERED = sorted(paired, key=lambda x: x[1])
            coord_lower_energy[ORDERED[0][0].split('.')[0]] = ORDERED[0][2]
            to_print = []
            for i in range(len(ORDERED)):
                to_print.append([ORDERED[i][0].split('_Cell_')[-1].split('.')[0], ORDERED[i][1]])
            to_print = pd.DataFrame(to_print, columns = ['Cell', 'Gibbs (hartree/partícula)'])    
            print(first_name)
            print(to_print)
            print('\n')
            with open(first_name+'.orden', 'wt') as final:
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
            index = [g.split('_Cell_')[-1].split('.')[0] for g in good]
            to_print_rmsd_matrix = pd.DataFrame(rmsd_matrix, index = index, columns = index)
            with open(first_name+'.rmsd', 'wt') as final:
                to_print_rmsd_matrix.to_string(final) 
        #========================Se exporta el .xyz==================
            if xyz_out:
                with open(ORDERED[0][0].split('.')[0]+'.xyz', 'wt') as final:
                    final.write(str(len(ORDERED[0][2]))+'\n\n')
                    ORDERED[0][2].to_string(final, header=False, index=False)

        
    #======================================================================
    if wrong_freq:
        print(f'Las cálculos: {wrong_freq} presentaron frecuencias negativas o imaginarias.')
    if wrong_end:
        print(f'Las cálculos: {wrong_end} presentaron problemas para un correcta finalizacion. Check the psi4 output!')

    if parameterize_path:
        for conf in coord_lower_energy:
            conf_path = os.path.join(parameterize_path, conf)
            tools.makedirs(conf_path)
            xyztmp = tempfile.NamedTemporaryFile(suffix='.xyz')
            with open(xyztmp.name, 'w') as t:
                t.write(f"{len(coord_lower_energy[conf])}\n\n")
                coord_lower_energy[conf].to_string(t, header=False, index=False)
            OBconvert.obconvert(xyztmp.name, os.path.join(conf_path, f"{conf}.mol2"))
            templates.PARAM(name = conf).write(os.path.join(conf_path, f"{conf}.sh"))

    return coord_lower_energy

 