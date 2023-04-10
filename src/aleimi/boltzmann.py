#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import rmsd
import numpy as np
import pandas as pd
import os
from aleimi import utils

def arc_reader (arc:str):
    """This , for sure, will change in the future to a more object-oriented paradigm
    It reads a arc MOPAC file and give a tuple of

    Parameters
    ----------
    arc : str
        The arc file path

    Returns
    -------
    list[tuple]
        A list of tuple sorted by energy. Each tuple contains:
            #. cells: conformer identifier,
            #. HeatsOfFormation_kcalmol: self-explanatory,
            #. CONTAINER__H: Coordinates numpy array of the heavy atoms with shape (Number of heavy atoms, 3),
            #. Class_E: The classic energy calculated during :meth:`aleimi.confgen.main`
    """
    with open(arc, 'rt', encoding='latin-1') as a:
        lines = a.readlines()

    # getting data from arc
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

def out_reader (out):
    """This , for sure, will change in the future to a more object-oriented paradigm
    It reads a out MOPAC file and give a tuple of

    Parameters
    ----------
    out : str
        The arc file path

    Returns
    -------
    list[tuple]
        A list of tuple sorted by energy. Each tuple contains:
            #. cells: conformer identifier,
            #. HeatsOfFormation_kcalmol: self-explanatory,
            #. CONTAINER__H: Coordinates numpy array of the heavy atoms with shape (Number of heavy atoms, 3)
            #. Class_E: The classic energy calculated during :meth:`aleimi.confgen.main`
    """
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
                    utils.ignoreLines(f, 1)
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

def main(file_path:str, Bd_rmsd:float = 1.0, Bd_E:float = 0.0, BOutPath:bool = True) -> pd.DataFrame:
    """It reads the MOPAC output file (arc or out) and create a Boltzmann table

    Parameters
    ----------
    file_path : str
        MOPAC output (.arc or .out). It will be read with :meth:`aleimi.boltzmann.arc_reader` or :meth:`aleimi.boltzmann.out_reader`
        depending on the extension.
    Bd_rmsd : float, optional
        RMSD to filter out redundant conformations, geometric filter, by default 1.0
    Bd_E : float, optional
        Energy difference in kJ to filter out redundant conformations, geometric filter, by default 1.0, by default 0.0
    BOutPath : bool, optional
        Directory to ouput the table {file_name}_boltzmann.csv", by default True

    Returns
    -------
    pd.DataFrame

    A Table with columns:
        #. `cell`: conformer identifier,
        #. `Class_E`: Classic energy from the RDKit optimization. in this example is `NaN` because we did not perform this optimization,
        #. `HeatOfFormation_kcal`: self-explanatory [kcal/mol],
        #. `Emin_Ei`: Difference in energy between the lower and the i-th conformer in [kcal/mol],
        #. `qi__Pi/Pmin__e^(Emin_Ei)/KbT`: Boltzmann factors,
        #. `Fraction_%__100*qi/q`: Occupancy of each conformer,

    Raises
    ------
    ValueError
        _description_
    """
    name, ext = os.path.splitext(os.path.basename(file_path))
    
    if ext == '.arc':
        ordered = (arc_reader(file_path))
    elif ext == '.out':
        ordered = (out_reader(file_path))
    else:
        raise ValueError(f"{file_path} does not have .arc or .out extension. Therefore is not readable by ALEIMI.")

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
    relative_kcalmol = [MinHeatsOfFormation_kcalmol - x for x in HeatsOfFormation_kcalmol]
    qi = [np.exp(E_r/(Kb*T)) for E_r in relative_kcalmol]
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
    DF['Emin_Ei'] = relative_kcalmol
    DF['qi__Pi/Pmin__e^(Emin_Ei)/KbT'] = qi
    DF['Fraction_%__100*qi/q'] = Fraction

    if BOutPath:
        with open(f"{name}_boltzmann.csv", 'wt') as rp:
            DF.to_csv(rp)
    
    return DF


if __name__ == '__main__':...