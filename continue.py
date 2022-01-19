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
from itertools import count
import os
from aleimi import templates
def main(directory, elapse):
    elapsed_paths = []
    elapsed_paths_tmp = []
    for root, dirs, files in os.walk(directory):
        if 'psi4.in' in files:
            if 'timer.dat' not in files:
                if len(elapsed_paths_tmp) <= elapse:
                    elapsed_paths_tmp.append(os.path.join(root, 'psi4.in'))
                else:
                    elapsed_paths.append(elapsed_paths_tmp)
                    elapsed_paths_tmp = []
    return elapsed_paths