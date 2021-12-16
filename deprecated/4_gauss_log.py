# -*- coding: utf-8 -*-
"""
===============================================================================
Created on    : Wed Sep  4 21:57:17 2019
Author        : Alejandro Martínez León
Mail          : [amleon@instec.cu, ale94mleon@gmail.com]
Affiliation   : Chemical Systems Modeling Group,
Affiliation   : Faculty of Radiochemistry, InSTEC-University of Havana, Cuba.
===============================================================================
DESCRIPTION   :
DEPENDENCIES  :
===============================================================================
"""
import glob
import rmsd
import numpy as np
import pandas as pd

# =============================================================================
#      USER SPECIFICATIONS
# =============================================================================

sh_files = False
sh_type_files = 'wahoo' # 'instec', 'wahoo'           

mem = '%mem=50000MB'
nproc = '%nprocshared=20'
keywords = '#P M062X/6-311++g(2df,2p) Integral=UltraFine nosymm density SCRF=(SMD, Solvent=Water) output=wfn scf=tight'
comments = 'Add comments here'
charge_mult = '0 1'
# =============================================================================


def no_repeat(l):
    result = []
    for item in l:
        if item not in result:
            result.append(item)
    return result

def log_read(log):
    
    global check_freq, Gibbs_free_E, exyz, xyz2RMSD_H
    with open(log, 'rt') as file:
        lines = file.readlines()
    check_freq = True
    exyz = []
    xyz2RMSD_H = []
    for i, line in enumerate(lines):
        if 'Frequencies' in line:
            freqs = line.split('--')[1].split()
            for freq in freqs:
                if float(freq) < 0:
                    check_freq = False
        if 'Sum of electronic and thermal Free Energies' in line:
            Gibbs_free_E = float(line.split('=')[1])
        if '!   Optimized Parameters   !' in line:
            stop = False
            for j in range(i, len(lines)):
                if 'orientation' in lines[j]:
                    for k in range(j+5, len(lines)):
                        if '-------' not in lines[k]:
                            exyz.append([int(lines[k].split()[1]), float(lines[k].split()[3]), float(lines[k].split()[4]), float(lines[k].split()[5])])
                            if int(lines[k].split()[1]) != 1:
                                xyz2RMSD_H.append([lines[k].split()[3], lines[k].split()[4], lines[k].split()[5]])   
                        else:
                            stop = True
                            break
                    if stop ==True:
                        break
    exyz = pd.DataFrame(exyz)
    xyz2RMSD_H = np.array(xyz2RMSD_H, dtype=np.float)
    return check_freq, Gibbs_free_E, exyz, xyz2RMSD_H 
        
        
        


if len(glob.glob('*.log')) == 0:
    exit('No se puede encontrar el archivo .log')

logs = glob.glob('*.log')
first_names = []
for log in logs:
    first_names.append(log.split('_Cell')[0])
first_names = no_repeat(first_names)


wrong_freq = []

for first_name in first_names:
    to_work = []
    for log in logs:
        if log.split('_Cell')[0] == first_name:
            to_work.append(log)           
    
    freq_checkes = []
    Gibbs_free_energies = []
    coords = []
    coord2RMSD_Hs =[]
    good = []
    for item in to_work:
        freq_check, Gibbs_free_energy, coord, coord2RMSD_H = (log_read(item))
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
                rmsd_matrix[i][j] = rmsd.kabsch_rmsd(P, Q, translate = True)
        
    to_print_rmsd_matrix = pd.DataFrame(rmsd_matrix, index = list(good), columns = list(good))
    with open(ORDERED[0][0].split('_Cell')[0]+'.rmsd', 'wt') as final:
        to_print_rmsd_matrix.to_string(final)
#========================Se exporta el .gjf==================
    with open(ORDERED[0][0].split('.')[0]+'_SP.gjf', 'wt') as final:
        final.write('%chk={}.chk\n'.format(ORDERED[0][0].split('.')[0]+'_SP'))
        final.write(mem+'\n')                        #
        final.write(nproc+'\n')                      #
        final.write(keywords+'\n\n')                 #
        final.write(comments+'\n\n')                 #
        final.write(charge_mult+'\n')                #
        ORDERED[0][2].to_string(final, header=False, index=False)
        final.write('\n \n')
        final.write(ORDERED[0][0].split('.')[0]+'_SP.wfn')
#======================================================================

#=====================se exporta el .sh==============================
# =============================================================================
#Obtencion de los archivos .sh
# =============================================================================
    if sh_files == True:
        if sh_type_files == 'wahoo': #'instec'   
            name = ORDERED[0][0].split('.')[0]+'_SP'
            with open(name+'.sh', 'wt') as final:
                final.write("#!/bin/bash\n\n")
                final.write("# bash a utiliser\n#$ -S /bin/bash\n\n")
                final.write("# Nom de la tache\n#$ -N %s\n\n" % (name))
                final.write("# Utiliser les variables d'environnement\n#$ -V\n\n" )     
                final.write("#$ -q normal.q\n\n")                      
                final.write("# je choisis l'environnement parallèle ( Parallel Environment ) (cf /drbd/sge/default/spool/qmaster/pe )\n#$ -pe shmem 12\n\n")                 
                final.write("# to launch g09 molecule gCD\ng09 < %s.gjf > %s.log\n\n" % (name,name))
        elif sh_type_files == 'instec':
            name = ORDERED[0][0].split('.')[0]+'_SP'
            numbthreads = int(nproc.split('=')[1])
            with open(name+'.sh', 'wt') as final:
                final.write("!/bin/bash\n\n")
                final.write("#$ -S /bin/bash\n\n")
                final.write("#$ -pe openmp %d\n\n" % (numbthreads))
                final.write("#$ -N %s\n\n" % (name))
                final.write("#$ -q all.q\n\n")
                final.write("#$ -cwd\n#$ -V\n\n" )
                final.write("#$ -m ea\n#$ -M amleon@instec.cu\n\n" )
                final.write("#$ -r y\n\n" )
                final.write("export OMP_NUM_THREADS=$NSLOTS\n\n")        
                final.write("cd $SGE_O_WORKDIR\n\n")   
                final.write("cp -r %s.gjf $TMPDIR\n\n" % (name))
                final.write("cd $TMPDIR\n\n")
                final.write("echo $OMP_NUM_THREADS >> screen.out\n")              
                final.write("g09 < %s.gjf > %s.log\n\n" % (name,name))
                final.write("cp -r * $SGE_O_WORKDIR\n\n")
            
#======================================================================
if len(wrong_freq) != 0:
    print('Las cálculos: %s presentaron frecuencias negativas.' % wrong_freq)
else:    
    print('Todas las frecuencias de los archivos en la carpeta de trabajo son mayores que cero.')
        
        
