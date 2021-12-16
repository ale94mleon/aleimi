#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
Created on    : Mon Dec 14 13:36:09 2020
Author        : Alejandro Martínez León
Mail          : [alejandro.marrtinezleon@uni-saarland.de, ale94mleon@gmail.com]
Affiliation   : Jochen Hub's Biophysics Group
Affiliation   : Faculty of NS, University of Saarland, Saarbrücken, Germany
===============================================================================
DESCRIPTION   :
DEPENDENCIES  :
===============================================================================
"""

import glob
import pandas as pd
from os import system

# =============================================================================
#      USER SPECIFICATIONS
#Si Energy es False generara tantos .gjf como se le señale en el variable
#molecules_cut, pero estas seran las mas favorables energeticamente.
# de desear extraerlas todas hacer molecules_cut = 0 y, por supesto, Energy = False
# Si se quiere obtener los archivos .sh tambien hacer True la variable sh_files. Sino
#hacerla False. en caso de desear modificar el archivo .sh ver el final de este codigo.
# =============================================================================
Energy = False
energy_cut = 1									
molecules_cut = 10                          
sh = True
mkdir = False #if true copy in subfolders       

level = 'MP2'
basis = '6-31+G(2d,2p)'
calc_type = 'VeryTightOpt'
freq = True
'''
This is a need for the Freq calcualtion in order to overcome probable memory shortage.
this use 8 Gb of memory per processor for the calculation. The user should provide the 
value according to the computer available memory. And only will be set if freq=True

The maximum memory on the smaug cluster  
'''
mem_per_mpi = 1200
SCF_details = 'VeryTightSCF'#S' Grid7 NoFinalGrid' #This is another need for freq calculation and will set only if freq = True
freq_type = 'NumFreq' #Could be AnFreq, or NumFreq, see manual of orca. only set if freq = True 
SMD = True
comments = '#! Computation at '+level+'/'+basis+' for '
charge_mult = '0 1'
#========================
partition = 'deflt'
mpi = 12
cpu_per_task = 1
time = '2-00:00'
mail_user = ''
number_nodes = 1
nice = 0
gpus = 0
exclude = ''#'fang[41-50]'
# =============================================================================

boltzmann_files = sorted(glob.glob('*.boltzmann'))
arc_files = sorted(glob.glob('*.arc'))

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

def toinsh(arc, to_extract, sh_file=False):
    with open(arc, 'rt', encoding='latin-1') as file:
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
            name = arc.split('.')[0]+'_Cell_'+str(cell)
            with open(name+'.inp', 'wt') as in_final:
                in_final.write(comments+name+'\n')
                in_final.write('%pal nprocs '+str(int(mpi))+' end\n')
                if freq: in_final.write('%maxcore '+str(mem_per_mpi)+'\n\n')
                in_final.write('! '+level+' '+basis+'\n')                        
                if freq: in_final.write('! '+SCF_details+'\n')
                in_final.write('! '+calc_type+'\n')
                if freq: in_final.write('! '+freq_type)
                if SMD: in_final.write(
'''
%cpcm 
	smd true
	SMDsolvent "water"
end
''')                   
                
                in_final.write('\n* xyz '+charge_mult+'\n')
                to_print.to_string(in_final, header=False, index=False)
                in_final.write('\n')
                in_final.write('*')
                print(name)

    # =============================================================================
    #Obtencion de los archivos .sh
    # =============================================================================
            if sh_file == True:
                with open(name+'.sh', 'wt') as sh_final:
                    sh_final.write('#!/bin/bash\n')
                    sh_final.write('#SBATCH --partition '+partition+' #which partition I want\n')
                    sh_final.write('#SBATCH --output '+name+'_myjob'+'.out #path for the slurm output\n')
                    sh_final.write('#SBATCH --error '+name+'_myjob'+'.err #path for the slurm error output\n')
                    sh_final.write('#SBATCH --mincpus='+str(mpi*cpu_per_task)+'\n')
                    sh_final.write('#SBATCH --ntasks='+str(mpi)+'\n')
                    sh_final.write('#SBATCH --cpus-per-task '+str(cpu_per_task)+' #number of cpu(logical cores)/task (task is normally an MPI process, default is one and the option to change it is -n)\n')
                    sh_final.write('##SBATCH --mem='+str(mpi*mem_per_mpi)+' #amount of real memory per allocated CPU required by the job.\n')
                    sh_final.write('#SBATCH --time '+time+' #how many time I want the resources (this impacts the job priority as well)\n')
                    sh_final.write('#SBATCH --job-name='+name+' #(to recognize your jobs when checking them with "squeue -u USERID")\n')
                    if mail_user: sh_final.write('#SBATCH --mail-user= '+mail_user+'\n')
                    sh_final.write('#SBATCH --nodes '+str(number_nodes)+' #number of node, usually 1 when no parallelization over nodes\n')
                    sh_final.write('#SBATCH --nice='+str(nice)+' #lowering your priority if >0\n')
                    sh_final.write('#SBATCH --gpus='+str(gpus)+' #number of gpu you want\n')
                    if exclude: sh_final.write('#SBATCH --exclude='+exclude+' #exclude some nodes\n')
                    sh_final.write(
'''
# This block is echoing some SLURM variables
echo "Jobid = $SLURM_JOBID"
echo "Host = $SLURM_JOB_NODELIST"
echo "Jobname = $SLURM_JOB_NAME"
echo "Subcwd = $SLURM_SUBMIT_DIR"
echo "SLURM_CPUS_PER_TASK = $SLURM_CPUS_PER_TASK"
                                    
# This block is for the execution of the program
export PATH="/data/shared/opt/ORCA/openmpi314/bin:$PATH"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/data/shared/opt/ORCA/openmpi314/lib"

# To order the folder
mkdir ${SLURM_JOB_NAME}
mv ${SLURM_JOB_NAME}* ${SLURM_JOB_NAME}
cd ${SLURM_JOB_NAME}

''')
                    sh_final.write('$(which orca) ${SLURM_JOB_NAME}.inp > ${SLURM_JOB_NAME}.log  #run command, this assumes that you change your .bash_profile to get orca') #--use-hwthread-cpus
# =============================================================================
#Programa principal. las moleculas a extraer se conforman de dos modos.
#print('Energy = %s\nenergy_cut = %d\nmolecules_cut = %d\nsh_files = %s\n' % (Energy, energy_cut, molecules_cut, sh_files))
# =============================================================================


to_extract = []
for i, boltzmann_file in enumerate(boltzmann_files):
    if Energy is True:
        cell = []
        df = pd.read_table(boltzmann_file, sep='\s+')
        df_subset = df[df.Emin_Ei >= -energy_cut]
        to_extract = df_subset.cell.tolist()    
        to_extract = sorted(to_extract)
        toinsh(arc_files[i],to_extract,sh_file=sh)
        
    elif molecules_cut != 0:
        cell = []
        df = pd.read_table(boltzmann_file, sep='\s+')
        df_subset = df.iloc[:molecules_cut]
        to_extract = df_subset.cell.tolist()    
        to_extract = sorted(to_extract)        
        toinsh(arc_files[i],to_extract,sh_file=sh)
    else:
        cell = []
        df = pd.read_table(boltzmann_file, sep='\s+')
        df_subset = df.iloc[:]
        to_extract = df_subset.cell.tolist()    
        to_extract = sorted(to_extract)
        
        toinsh(arc_files[i],to_extract,sh_file=sh)


if mkdir:
    names = [item.split('.')[0] for item in glob.glob('*.in')]
    for name in names:
        system('mkdir '+name)
        system('mv '+name+'.* '+name)
