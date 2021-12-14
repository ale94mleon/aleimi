#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 11:21:11 2019

@author: alejandro
"""
#toma una trayectoria y la transforma en .mop
import pandas as pd
import glob as glob


if len(glob.glob('*.HCS')) == 0:
    exit('No se puede encontrar el archivo .HCS')
# =============================================================================
#      USER SPECIFICATIONS
# =============================================================================
#cut = 308                            
hcss = glob.glob('*.HCS')

keywords = 'PM7 precise ef xyz geo-ok t=3d'
# =============================================================================


for hcs in hcss:
    

    cont=0
    
    with open(hcs, 'rt') as file:
        lines = file.readlines()
        
    final=open(hcs.split('.')[0]+'.mop','w')
    
    atoms=[]
    for i, line in enumerate(lines):
        if 'endmol 1' in line:
            natoms = int(lines[i-1].split()[1])
            chunk1 = lines[(i-(natoms)):(i)]
            for c in chunk1:
                atoms.append(c.strip().split()[3])
        if ('[Conformation' in line):#and cont<cut:
            cont+=1
            cell = int(line.strip().split()[-1][:-1])
            Class_E = float(lines[i+1].split('=')[1])
            chunk=lines[i+7:i+7+natoms]
            X = []
            Y = []
            Z =[]
            for c in chunk:
                X.append(float((c.strip().split('=')[1].split()[0])))
                Y.append(float((c.strip().split('=')[1].split()[1])))
                Z.append(float((c.strip().split('=')[1].split()[2])))
            d = {'Element':atoms, 'X':X, 'Y':Y, 'Z':Z}
            df = pd.DataFrame(data=d)    
            to_print = df[['Element', 'X', 'Y', 'Z']]
            final.write(keywords+'\n')                 
            final.write('E = ' +str(Class_E)+' Molécula en el opt: '+ str(cell)+'\n')                 
            final.write('CELL: %d\n' % (cell))
            to_print.to_string(final, header=False, index=False)
            final.write('\n0\n')
