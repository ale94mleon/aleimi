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
energy_cut = 2									# Reading from .boltzmann (kcal/mol)
molecules_cut = 1                          
sh = True
mkdir = True           

mem = '10 GB'
level = 'b3lyp'
basis = '6-31+G(d,p)'
calc_type = 'optimize'
freq = True
PCM = False
comments = '#! Computation at '+level+'/'+basis+' for '
charge_mult = '0 1'
#========================
partition = 'deflt'
nproc = 12
time = '2-00:00'
mail_user = ''
number_nodes = 1
nice = 0
gpus = 0
exclude = 'fang[41-50]'
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
            with open(name+'.in', 'wt') as in_final:
                in_final.write(comments+name+'\n\n')
                in_final.write('memory '+mem+'\n\n')                        
                in_final.write('molecule '+name+' {\n')
                in_final.write(charge_mult+'\n')
                to_print.to_string(in_final, header=False, index=False)
                in_final.write('\n')
                if PCM: in_final.write('symmetry c1\n')
                in_final.write('}\n\n')
                if PCM:
                    in_final.write('set {\n  basis '+basis)
                    in_final.write(
'''
  scf_type pk
  pcm true
  pcm_scf_type total
}

pcm = {
   Units = Angstrom
   Medium {
   SolverType = IEFPCM
   Solvent = Water
   }

   Cavity {
   RadiiSet = UFF
   Type = GePol
   Scaling = False
   Area = 0.3
   Mode = Implicit
   }
}
''')
                else:                
                    in_final.write('set basis '+basis+'\n')
                in_final.write(calc_type+"(\'"+level+"\')\n")
                if freq: in_final.write(level+"_e, "+level+"_wfn = frequencies"+"(\'"+level+"\', return_wfn=True, dertype=1)")
                print(name)

    # =============================================================================
    #Obtencion de los archivos .sh
    # =============================================================================
            if sh_file == True:
                with open(name+'.sh', 'wt') as sh_final:
                    sh_final.write('#!/bin/bash\n')
                    sh_final.write('#SBATCH -p '+partition+' #which partition I want\n')
                    sh_final.write('#SBATCH -o myjob_'+name+'.out #path for the slurm output\n')
                    sh_final.write('#SBATCH -e myjob_'+name+'.err #path for the slurm error output\n')
                    sh_final.write('#SBATCH -c '+str(nproc)+' #number of cpu(logical cores)/task (task is normally an MPI process, default is one and the option to change it is -n)\n')
                    sh_final.write('#SBATCH -t '+time+' #how many time I want the resources (this impacts the job priority as well)\n')
                    sh_final.write('#SBATCH --job-name='+name+' #(to recognize your jobs when checking them with "squeue -u USERID")\n')
                    sh_final.write('#SBATCH --mail-user= '+mail_user+'\n')
                    sh_final.write('#SBATCH -N '+str(number_nodes)+' #number of node, usually 1 when no parallelization over nodes\n')
                    sh_final.write('#SBATCH --nice='+str(nice)+' #lowering your priority if >0\n')
                    sh_final.write('#SBATCH --gpus='+str(gpus)+' #number of gpu you want\n')
                    sh_final.write('#SBATCH --exclude='+exclude+' #exclude some nodes\n')
                    sh_final.write(
'''
# This block is echoing some SLURM variables
echo "Jobid = $SLURM_JOBID"
echo "Host = $SLURM_JOB_NODELIST"
echo "Jobname = $SLURM_JOB_NAME"
echo "Subcwd = $SLURM_SUBMIT_DIR"
echo "SLURM_CPUS_PER_TASK = $SLURM_CPUS_PER_TASK"
                                    
# This block is for the execution of the program
source /home/users/all-jh/opt/miniconda3/etc/profile.d/conda.sh #It is needed to use the conda activate command
conda activate htmd
#In order to have a temp file with more capacity than /temp (default)
cd $SLURM_SUBMIT_DIR
mkdir psi4.$SLURM_JOBID
MYSCRATCH=$SLURM_SUBMIT_DIR/psi4.$SLURM_JOBID
export PSI_SCRATCH=$MYSCRATCH\n
''')
                    sh_final.write('psi4 -i ${SLURM_JOB_NAME}.in -o ${SLURM_JOB_NAME}.out -n '+str(nproc)+' #run command, in, out and number of threads to be used\n\n')
                    sh_final.write('rm -rf $MYSCRATCH')
# =============================================================================
#Programa principal. las moleculas a extraer se conforman de dos modos.
#print('Energy = %s\nenergy_cut = %d\nmolecules_cut = %d\nsh_files = %s\n' % (Energy, energy_cut, molecules_cut, sh_files))
# =============================================================================
print('Las variables predefinidas fueron:')
if Energy==True and sh==True:
    print('energy_cut = %d\nsh_files = %s\nsh_type_files = slrum_psi4\n' % (energy_cut, sh))
    print('energy_cut = %d\n' % (energy_cut))
elif Energy==False and sh==True:
    print('molecules_cut = %d\nsh_files = %s\nsh_type_files = slrum_psi4\n' % (molecules_cut, sh))
elif Energy==False and sh==False:
    print('molecules_cut = %d\n' % (molecules_cut))




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
