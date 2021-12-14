# -*- coding: utf-8 -*-
"""
===============================================================================
Created on    : Tue Aug 27 14:05:37 2019
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
import pandas as pd

# =============================================================================
#      USER SPECIFICATIONS
#Si Energy es False generara tantos .gjf como se le señale en el variable
#molecules_cut, pero estas seran las mas favorables energeticamente.
# de desear extraerlas todas hacer molecules_cut = 0 y, por supesto, Energy = False
# Si se quiere obtener los archivos .sh tambien hacer True la variable sh_files. Sino
#hacerla False. en caso de desear modificar el archivo .sh ver el final de este codigo.
# =============================================================================
Energy = False
energy_cut = 2									# Reading from .boltzmann (kcal/mol)
molecules_cut = 5                           
sh_files = True
sh_type_files = 'wahoo' # 'instec', 'wahoo'           

mem = '%mem=50000MB'
nproc = '%nprocshared=20'
keywords = '# opt freq M062X/6-31+g(2d,2p) Integral=UltraFine nosymm SCRF=(SMD, Solvent=Water)'
comments = 'Add comments here'
charge_mult = '0 1'
# =============================================================================

boltzmann_files = glob.glob('*.boltzmann')
arc_files = glob.glob('*.arc')

# =============================================================================
#     Comprobacion de que los archivos .arc y .boltzmann son correctos
# =============================================================================


test_boltzmann = []
for item in boltzmann_files:
    test_boltzmann.append(item.split('.')[0])

test_arc = []
for item in arc_files:  
    test_arc.append(item.split('.')[0])

if test_boltzmann==test_arc:
    print('Los archivos .arc y .boltzmann son correctos\n')
else:
    exit('Problemas con los archivos .arc y/o .boltzmann. No hay la misma cantidad o tienen diferentes nombres.')

# =============================================================================
#Toma un archivo .arc y extrae las celdas que esten dentro de to_extract
#y estas celdas las convierte en .gjf
# =============================================================================

def togjf(arc, to_extract):
    with open(arc, 'rt') as file:
        lines = file.readlines()
    
    for line in lines:
        if 'Empirical Formula' in line:
            natoms = int(line.split()[-2])
            break
    for i, line in enumerate(lines):
        if ('FINAL GEOMETRY OBTAINED' in line) and (int(lines[i+3].strip().split(':')[1]) in to_extract):
            cell = int(lines[i+3].strip().split(':')[1])
            chunk = lines[i+4:i+4+natoms]
            sliced = []
            for c in chunk:
                sliced.append(c.split())
            df = pd.DataFrame(sliced)
            to_print = df[[0, 1, 3, 5]]
    
            with open(arc.split('.')[0]+'_Cell_'+str(cell)+'.gjf', 'wt') as final:
                final.write('%chk={}.chk\n'.format(arc.split('.')[0]+'_Cell_'+str(cell)))
                final.write(mem+'\n')                        #
                final.write(nproc+'\n')                      #
                final.write(keywords+'\n\n')                 #
                final.write(comments+'\n\n')                 #
                final.write(charge_mult+'\n')                #
                to_print.to_string(final, header=False, index=False)
                final.write('\n ')
                print(arc.split('.')[0]+'_Cell_'+str(cell))


# =============================================================================
#Programa principal. las moleculas a extraer se conforman de dos modos.
#print('Energy = %s\nenergy_cut = %d\nmolecules_cut = %d\nsh_files = %s\n' % (Energy, energy_cut, molecules_cut, sh_files))
# =============================================================================
print('Las variables predefinidas fueron:')
if Energy==True and sh_files==True:
    print('energy_cut = %d\nsh_files = %s\nsh_type_files = %s\n' % (energy_cut, sh_files, sh_type_files))
elif Energy==True and sh_files==False:
    print('energy_cut = %d\n' % (energy_cut))
elif Energy==False and sh_files==True:
    print('molecules_cut = %d\nsh_files = %s\nsh_type_files = %s\n' % (molecules_cut, sh_files, sh_type_files))
elif Energy==False and sh_files==False:
    print('molecules_cut = %d\n' % (molecules_cut))



print('Se extrajeron las celdas:')

to_extract = []
for i, boltzmann_file in enumerate(boltzmann_files):
    if Energy is True:
        cell = []
        df = pd.read_table(boltzmann_file, sep='\s+')
        df_subset = df[df.Emin_Ei >= -energy_cut]
        to_extract = df_subset.cell.tolist()    
        to_extract = sorted(to_extract)
        togjf(arc_files[i],to_extract)
        
    elif molecules_cut != 0:
        cell = []
        df = pd.read_table(boltzmann_file, sep='\s+')
        df_subset = df.iloc[:molecules_cut]
        to_extract = df_subset.cell.tolist()    
        to_extract = sorted(to_extract)
        
        togjf(arc_files[i],to_extract)
    else:
        cell = []
        df = pd.read_table(boltzmann_file, sep='\s+')
        df_subset = df.iloc[:]
        to_extract = df_subset.cell.tolist()    
        to_extract = sorted(to_extract)
        
        togjf(arc_files[i],to_extract)
        

        
# =============================================================================
#Obtencion de los archivos .sh
# =============================================================================
if sh_files == True:
    if sh_type_files == 'wahoo': #'instec'   
        fields=glob.glob('*.gjf')
        for item in fields:
            name = item.split('.')[0]
            with open(name+'.sh', 'wt') as final:
                final.write("#!/bin/bash\n\n")
                final.write("# bash a utiliser\n#$ -S /bin/bash\n\n")
                final.write("# Nom de la tache\n#$ -N %s\n\n" % (name))
                final.write("# Utiliser les variables d'environnement\n#$ -V\n\n" )     
                final.write("#$ -q normal.q\n\n")                      
                final.write("# je choisis l'environnement parallèle ( Parallel Environment ) (cf /drbd/sge/default/spool/qmaster/pe )\n#$ -pe shmem 12\n\n")                 
                final.write("# to launch g09 molecule gCD\ng09 < %s.gjf > %s.log\n\n" % (name,name))
    elif sh_type_files == 'instec':
        fields=glob.glob('*.gjf')
        for item in fields:
            name = item.split('.')[0]
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