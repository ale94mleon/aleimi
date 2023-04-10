# -*- coding: utf-8 -*-
"""
===============================================================================
Created on    : Thu Aug 22 22:23:48 2019
Author        : Alejandro Martínez León
Mail          : [amleon@instec.cu, ale94mleon@gmail.com]
Affiliation   : Chemical Systems Modeling Group,
Affiliation   : Faculty of Radiochemistry, InSTEC-University of Havana, Cuba.
===============================================================================
DESCRIPTION   :
DEPENDENCIES  :
===============================================================================
"""

import rmsd
import numpy as np
import pandas as pd
import os

def arc_reader (arcf):
    with open(arcf, 'rt', encoding='latin-1') as a:
        lines = a.readlines()

    # getting data from arc                                                       #
    HeatsOfFormation_kcalmol = []
    cells = []
    Class_E = []
    # finding No. of atoms
    for line in lines:
        if 'Empirical Formula' in line:
            natoms = int(line.split()[-2])
            break
    
    #CONTAINER = []
    CONTAINER__H = []
    #cart = []
    cart__H = []
    #atoms = []
    for i, line in enumerate(lines):

        if 'HEAT OF FORMATION' in line:
            #I am getting the one in kcal/mol
            HeatsOfFormation_kcalmol.append(float(line.split()[-5])) #HEAT OF FORMATION       =       -180.69875 KCAL/MOL =    -756.04358 KJ/MOL
       
        elif 'Empirical Formula' in line:
            try:
                Class_E.append(float(lines[i+3].split('=')[-1].split()[0])) #las dos versiones con y sin optimizacion
            except:
                Class_E.append(None)
                    
            cells.append(int(lines[i+4].split(':')[1]))
    
        elif ('FINAL GEOMETRY OBTAINED' in line):
            chunk = lines[i+4:i+4+natoms]
            cart__H = []
            for c in chunk:
                #atoms.append(c.split()[0])
                #cart.append([c.split()[0].strip(), float(c.split()[1]), float(c.split()[3]), float(c.split()[5])])
                if c.split()[0] != 'H':
                    cart__H.append([c.split()[1], c.split()[3], c.split()[5]]) # No estoy tomando los atomos, solamente coordenadas c.split()[0].strip(),
            #CONTAINER.append(np.asarray(pd.DataFrame(cart)))
            CONTAINER__H.append(np.array(cart__H, dtype=np.float64))
    #atoms = (np.asarray(atoms))
    # .... organizing
    paired = list(zip(cells, HeatsOfFormation_kcalmol, CONTAINER__H, Class_E)) # Esto genera un arreglo de tuplas, me une los arreglos
    ORDERED = sorted(paired, key=lambda x: x[1])  #Esto ordena la tupla segun la energia de menor a mayor
    return ORDERED #, atoms]

def ignoreLines(f, n):
    for _ in range(n): f.readline()

def out_reader (out):
    f = open(out, 'r')
    chunk = []
    HeatsOfFormation_kcalmol = []
    cells = []
    Class_E = []
    CONTAINER__H = []
    cart__H = []


    # getting data from out                                                       #
    # finding No. of atoms
    while True:
        line = f.readline()
        if "Empirical Formula" in line:
            natoms = int(line.split()[-2])
            break

    f = open(out, 'r')
    while True:
        line = f.readline()
        if len(line) == 0:break
               
        if 79*'-' in line:
            while True:
                line = f.readline()
                if (79*'*' in line) or (len(line) == 0):break
                if "E_UFF = " in line:
                    split = line.split('=')[-1]
                    try:
                        Class_E.append(float(split))
                    except:
                        Class_E.append(None)

                elif 'CELL' in line:
                    cells.append(int(line.split(':')[1]))
                    
                elif 'HEAT OF FORMATION' in line:
                    HeatsOfFormation_kcalmol.append(float(line.split()[-5]))
                elif 'CARTESIAN COORDINATES' in line:
                    ignoreLines(f, 1)
                    cont = 0
                    chunk = []
                    while cont < natoms:
                        chunk.append(f.readline())
                        cont += 1
                    cart__H = []
                    for c in chunk:
                        if c.split()[1] != 'H':
                            cart__H.append([c.split()[2], c.split()[3], c.split()[4]]) # No estoy tomando los atomos, solamente coordenadas c.split()[0].strip(),

            CONTAINER__H.append(np.array(cart__H, dtype=np.float64))
    f.close
    # .... organizing
    paired = list(zip(cells, HeatsOfFormation_kcalmol, CONTAINER__H, Class_E)) # Esto genera un arreglo de tuplas, me une los arreglos
    ORDERED = sorted(paired, key=lambda x: x[1])  #Esto ordena la tupla segun la energia de menor a mayor
    return ORDERED

def main(file_path, Bd_rmsd = 1.0, Bd_E = 0.0, BOutPath = True):
    """This function generates a table with the population of each conformer respect to a Boltzmann probability distribution.

    Args:
        file_path (path): arc or out file. It will be read with arc_reader or out_reader depending on the extension.
        Bd_rmsd (float): Difference in RMSD (Angstrom) between the conformers. Defaults is 1.0.
        Bd_E (float, optional): In case that also you need to screen with a threshold of energy, you specify the value in kcal/mol. Defaults to None.
        BOutPath (path): In case that you want to write down the table, the path for the file.

    Raises:
        ValueError: The Boltzmann table
    """

    # d_E = 0.001
    # d_rmsd = 1.0
    ext = os.path.basename(file_path).split('.')[-1]
    name = os.path.basename(file_path)[:-(len(ext) + 1)]
    
    if ext == 'arc':
        ordered = (arc_reader(file_path))
    elif ext == 'out':
        ordered = (out_reader(file_path))
    else:
        raise ValueError(f"{file_path} does not have .arc or .out extension. Therefore is not readeable by ALEIMI.")

    if Bd_E:
        for i, _ in enumerate(ordered):
            to_trash_degenerated = []
            for idx, _ in enumerate(ordered):
                if i < idx:
                    Ei = ordered[i][1]
                    Eidx = ordered[idx][1]
                    delta = abs(Ei - Eidx)
                    if delta <= Bd_E:
    # =============================================================
    #     CHECKING Geometric degeneracy
    # =============================================================
                        P = ordered[i][2]
                        Q = ordered[idx][2]

                        RMSD = rmsd.kabsch_rmsd(P, Q, translate = True)

                        if RMSD <= Bd_rmsd:
                        # reject identical structure
                            to_trash_degenerated.append(idx)
# =========================================================================
#     FOR EACH STRUCTURE, eliminate degenerated and save lot of time
# =========================================================================
            to_trash_degenerated = sorted(to_trash_degenerated, reverse=True)
            _ = [ordered.pop(x) for x in to_trash_degenerated]

    else:
        for i, x in enumerate(range(len(ordered))):
            to_trash_degenerated = []
            for idx, y in enumerate(range(len(ordered))):
                if i < idx:
                
# =============================================================
#     CHECKING Geometric degeneracy
# =============================================================
                    P = ordered[i][2]
                    Q = ordered[idx][2]

                    RMSD = rmsd.kabsch_rmsd(P, Q, translate = True)

                    if RMSD <= Bd_rmsd:
                        # reject identical structure and kept the lowest energy (because ordered() is ordered using the enrgy, so idx always will have a grater enrgy)
                        to_trash_degenerated.append(idx)

# =========================================================================
#     FOR EACH STRUCTURE, eliminate degenerated and save lot of time
# =========================================================================
            to_trash_degenerated = sorted(to_trash_degenerated, reverse=True)
            [ordered.pop(x) for x in to_trash_degenerated]
        


# =============================================================================
#      WORKING with UNDEGENERATED. Cambie la manera de calculos los parametros:
#Me base en: James B. Foresman - Exploring Chemistry With Electronic Structure Methods 3rd edition (2015) pag 182
# y Mortimer_Physical Chemistry_(3rd.ed.-2008) pag 1045
# =============================================================================
    Kb = 1.987204259E-3                        # kcal/(mol⋅K)
    T = 298.15                                 # Absolute T (K)
    DF = pd.DataFrame()
    cells = [ordered[i][0] for i, x in enumerate(ordered)]
    Class_E = [ordered[i][3] for i, x in enumerate(ordered)]
    HeatsOfFormation_kcalmol = [ordered[i][1] for i, x in enumerate(ordered)]
    MinHeatsOfFormation_kcalmol = min(HeatsOfFormation_kcalmol)
    relative_kJ = [MinHeatsOfFormation_kcalmol - x for x in HeatsOfFormation_kcalmol]
    qi = [np.exp(E_r/(Kb*T)) for E_r in relative_kJ]
    q = sum(qi)
    Fraction = [100*i/q for i in qi]
    #Z = [np.e**(-(E/(k*T))) for E in energy_kcal] #no pudo calcular Z: verflowError: (34, 'Result too large')
    #Pi_b = [(np.e**-(E/(k*T)))/Z for E in energy_kcal]
    # =============================================================================
    #     DATAFRAME
    # =============================================================================
    DF['cell'] = cells
    DF['Class_E'] = Class_E
    DF['HeatOfFormation_kcal/mol'] = HeatsOfFormation_kcalmol
    DF['Emin_Ei'] = relative_kJ
    DF['qi__Pi/Pmin__e^(Emin_Ei)/KbT'] = qi
    DF['Fraction_%__100*qi/q'] = Fraction

    if BOutPath:
        with open(f"{name}_boltzmann.csv", 'wt') as rp:
            DF.to_csv(rp)
    return DF



