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
import os
from aleimi import templates, utils
def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def main(directory, elapse, launch = False):
    all_paths = []
    for root, dirs, files in os.walk(directory):
        root = os.path.abspath(root)
        if 'psi4.in' in files:
            if 'timer.dat' not in files:
                all_paths.append(root)

    for (i,elapsed_paths) in enumerate(chunks(all_paths, elapse)):
        T = templates.CONTINUE(elapsed_paths, machine = 'gwdg', name = f"elapse{i+1}")
        T.write(f"elapse{i+1}.sh")
        if launch == True:
            utils.run(f"sbatch elapse{i+1}.sh")
    

    
    
