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
import os


def ignoreLines(f, n):
    for i in range(n): f.readline()

#def extract(arc_file, boltzmann_file, energy_cut = 2, conformer_cut = None, mksh = True, mkdir = True):
# Energy cut is in kcal/mol
def extract(boltzmann_file, energy_cut = 2, conformer_cut = None):
    if energy_cut:
        df = pd.read_table(boltzmann_file, sep='\s+')
        df_subset = df[df.Emin_Ei >= -energy_cut]
        indx_to_extract = df_subset.cell.tolist()    
        indx_to_extract = sorted(indx_to_extract)

        
    elif conformer_cut:
        df = pd.read_table(boltzmann_file, sep='\s+')
        df_subset = df.iloc[:conformer_cut]
        indx_to_extract = df_subset.cell.tolist()    
        indx_to_extract = sorted(indx_to_extract)        

    else:
        df = pd.read_table(boltzmann_file, sep='\s+')
        df_subset = df.iloc[:]
        indx_to_extract = df_subset.cell.tolist()    
        indx_to_extract = sorted(indx_to_extract)        
    return indx_to_extract

def get_coords(file_path, indx_to_extract):
    file_name = os.path.basename(file_path).split('.')[0]
    file_ext = os.path.basename(file_path).split('.')[-1]

    if file_ext == 'arc':
        with open(file_path, 'rt', encoding='latin-1') as file:
            lines = file.readlines()
        
        to_return = []
        for line in lines:
            if 'Empirical Formula' in line:
                natoms = int(line.split()[-2])
                break
        for i, line in enumerate(lines):
            if ('FINAL GEOMETRY OBTAINED' in line) and (int(lines[i+3].strip().split(':')[1]) in indx_to_extract):
                cell = int(lines[i+3].strip().split(':')[1])
                chunk = lines[i+4:i+4+natoms]
                sliced = []
                for c in chunk:
                    sliced.append(c.split())
                df = pd.DataFrame(sliced)
                # Creating a tuple (conf_name_Cell, coords)
                to_return.append((f"{file_name}_Cell_{cell}", df[[0, 1, 3, 5]]))
        return to_return

    elif file_ext == '.out':
        f = open(file_path, 'r')
        chunk = []
        cart = []


        # getting data from out                                                       #
        # finding No. of atoms
        
        while True:
            l = f.readline()
            if "Empirical Formula" in l:
                natoms = int(l.split()[-2])
                f.close()
                break
                            
        f = open(file_path, 'r')
        to_return = []
        while True:
            l = f.readline()
            k = l
            if len(l) == 0:break
                
            if 79*'-' in l:
                while True:
                    k = f.readline()
                    if (79*'*' in k) or (len(k) == 0):break

                    if 'CELL' in k: 
                        cell = int(k.split(':')[1])

                    elif 'CARTESIAN COORDINATES' in k and cell in indx_to_extract:
                        ignoreLines(f, 1)
                        cont = 0
                        chunk = []        
                        while cont < natoms:
                            chunk.append(f.readline())
                            cont += 1
                        cart = []
                        for c in chunk:
                            cart.append(c.split())
                            df = pd.DataFrame(cart)
                            to_return.append((f"{file_name}_Cell_{cell}",df[[1, 2, 3, 4]]))  
        return to_return
    else:
        raise ValueError(f"{file_path} does not have .arc or .out extension. Therefore is not readeable by ALEIMI.")

