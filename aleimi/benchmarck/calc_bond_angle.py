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

import numpy as np
import math as mt
from rmsd import kabsch_rmsd
import math as mt
def calc_dihedral(p1, p2, p3, p4):
    """ Calculate dihedral angle method.
    The angle is in [-pi, pi].
    """
    #p1 = np.array(p1)
    #p2 = np.array(p2)
    #p3 = np.array(p3)
    #p4 = np.array(p4)

    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3
    
    n1 = np.cross(b1,b2)
    n1 = n1/np.linalg.norm(n1)
    
    n2 = np.cross(b2,b3)
    n2 = n2/np.linalg.norm(n2)
    
    b2 = b2/np.linalg.norm(b2)
    
    m1 = np.cross(n1,b2)
    
    x = (n2*n1).sum()
    y = (n2*m1).sum()
    
    return mt.atan2(y,x)

def calc_angle_bond(p1, p2, p3):
    """ Calculate bond angle method.
    The angle is in [-pi, pi].
    """
    #p1 = np.array(p1)
    #p2 = np.array(p2)
    #p3 = np.array(p3)    
    b1 = p1 - p2
    b2 = p3 - p2
    b1b2 = (b1*b2).sum()
    norm_b1 = np.linalg.norm(b1)
    norm_b2 = np.linalg.norm(b2)
    return mt.acos(b1b2/(norm_b1*norm_b2))

def calc_bond_dist(p1, p2):
    #p1 = np.array(p1)
    #p2 = np.array(p2)    
    return np.linalg.norm(p1 - p2)

def psi4out2coords(out, element = False):
    '''
    

    Parameters
    ----------
    out : string
        name of the psi4.out file from opt.
    element : TYPE, optional
        DESCRIPTION. The default is False.
        if False:
            return only the xyz coords
        if True:
            Also return the element
    Returns
    -------
    None.

    '''
    with open(out, 'rt', encoding='latin-1') as file:
        lines = file.readlines()
   
    xyz = []

    for i in range(len(lines)):
        if 'Final optimized geometry and variables:' in lines[i]:
            for j in range(i+6,len(lines)):
                if '\n' == lines[j]: break
                split = lines[j].split()    
                if element:
                    xyz.append(np.array([split[0], float(split[1]), float(split[2]), float(split[3])]))
                else:
                    xyz.append(np.array([float(split[1]), float(split[2]), float(split[3])]))
    return xyz

#print(psi4out2coords('cis1_bh267_m_vacumm_Cell_682.out'))

def RMSD_xyz(out1, out2, H = False):
    '''
    

    Parameters
    ----------
    out1 : str
        name of the first psi4out.
    out2 : str
        ame of the second psi4out.
    H : TYPE, bool
        if true use the H atoms for the calculation. The default is False.

    Returns
    -------
    float
        rmsd_kabsch with translate = True

    '''
    exyz1 = psi4out2coords(out1, element = True)
    exyz2 = psi4out2coords(out2, element = True)
    if H:
        exyz1 = np.array([item[1:] for item in exyz1], dtype=np.float)
        exyz2 = np.array([item[1:] for item in exyz2], dtype=np.float)
    else:    
        exyz1 = np.array([item[1:] for item in exyz1 if item[0] != 'H'], dtype=np.float)
        exyz2 = np.array([item[1:] for item in exyz2 if item[0] != 'H'], dtype=np.float)
    return float(kabsch_rmsd(exyz1, exyz2, translate = True))

def psi4_freq_check(out):
    '''
    

    Parameters
    ----------
    out : str
        name of the psi4 out file.

    Returns
    -------
    bool
        True if the freq calculation was made and any imaginary modes are present.
        False if not
    '''
    with open(out, 'rt', encoding='latin-1') as file:
        lines = file.readlines()
    cont = 0 #To check if the freq calcualtion was made
    for line in lines[::-1]: #The freq calculation is at the end, less line to process
        if 'Freq [cm^-1]' in line:
            cont += 1
            freqs = line.split()[2:]
            for f in freqs:
                '''
                 psi4 represent the imaginary number as 5i, 
                 if an error ocurreduring the float convertion, 
                 or it was converted but if less than cero (another error on sqrt) 
                 then return False
                '''
                try:
                    mt.sqrt(float(f))
                except:
                    return False
    if cont > 0:
        return True
    else:
        return False

def psi4_opt_check(out):
    '''
    

    Parameters
    ----------
    out : str
        name of the psi4 out file.

    Returns
    -------
    bool
        True if the opt was a succes.

    '''
    with open(out, 'rt', encoding='latin-1') as file:
        lines = file.readlines()
    for line in lines:
        if 'Optimization is complete!' in line: return True
    return False



def transpouse(array):
    '''
    

    Parameters
    ----------
    array : a list of list
        DESCRIPTION.

    Returns
    -------
    inv : numpy array
        the first row now is the first columns
        the second row now is the second columns
        etc.

    '''
    array = np.asarray(array)
    inv = np.zeros(len(array)*len(array[0]), dtype=array.dtype)
    inv.shape = (len(array[0]), len(array))
    for i in range(len(array)):
        for j in range(len(array[0])):
            inv[j,i]=array[i,j]
    return inv

def errors(dif):
    '''
    

    Parameters
    ----------
    dif : list of numbers (array)
        is a list that contain the diferences between the test set of values
        and the reference.
        It is crucial the order because the interpretation change
        test - reference is used in this function

    Returns
    -------
    RMSE : float
        Root Mean Square Error.
    MSE : float
        Mean Signed Error.
    MUE : float
        Mean Unsigned Error.

    '''
    RMSE = mt.sqrt(sum([item**2 for item in dif])/len(dif))
    MSE = sum(dif)/len(dif)
    MUE = sum([abs(item) for item in dif])/len(dif)
    return [RMSE, MSE, MUE]

def psi4_free_energy(out):
    with open(out, 'rt', encoding='latin-1') as file:
        lines = file.readlines()
    for line in lines[::-1]: #a[start:stop:step]
        if 'Total G, Free enthalpy at  298.15 [K]' in line:
            return float(line.split('[K]')[1].split()[0])
