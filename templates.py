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
import copy, os
class INPUT:
    def __init__(self, engine, machine = 'smaug', **keywords):
        self.engine = engine
        self.machine = machine
        self.default_partition = {'smaug': 'deflt', 'gwdg': 'medium'}
        self.default_exclude = {'smaug': "fang[1,11-50]", 'gwdg': None}
        self.default_ntasks = {'smaug': 12, 'gwdg': 24}
        self.default_keywords = {
            'psi4':{
                #Parameters for input
                'theory': 'MP2',
                'basis': '6-31+G*',
                'memory': '10 GB',
                'name': 'mol',
                'charge': 0,
                'multiplicity': 1,
                'PCM': False,
                'coords': pd.DataFrame(),
                'calc_type': 'optimize',
                'freq': True,
                #Parameters for jobsh
                'time': '2-00:00',
                'ntasks': self.default_ntasks[self.machine],
                'cpus-per-task': 1,
                'partition': self.default_partition[self.machine],
                'nodes': 1,
                'nice': 0,
                'gpus': 0,
                'mail_user': None,
                'exclude': self.default_exclude[self.machine],
                },
            'orca':{
                'theory': 'RI-MP2',
                'basis': '6-31+G* def2-TZVP/C',
                'name': 'mol',
                'charge': 0,
                'multiplicity': 1,
                'SMD': False,
                'coords': pd.DataFrame(),
                'calc_type': 'Opt',
                'freq_type': 'NumFreq',
                'SCF_details': 'VeryTightSCF VeryTightOpt',
                'maxcore':2400,
                #Parameters for jobsh
                'time': '2-00:00',
                'ntasks': 6,
                'cpus-per-task': 2,
                'partition': 'deflt',
                'nodes': 1,
                'nice': 0,
                'gpus': 0,
                'mail_user': None,
                'exclude': "fang[1,11-50]" #None
                },
            'gaussian':{
                'theory': 7,
                'basis': 2,
                'memory':4,
                'name':5,
                'charge':0,
                'multiplicity':1,
                'PCM': False,
                'coords': None,
                'calc_type':5,
                'freq': False}
        }
        self.keywords = copy.deepcopy(self.default_keywords)
        for key in keywords:
            self.keywords[engine][key] = keywords[key]
        self.input = ''
        self.jobsh = ''
        self.parse()
    
    def parse(self):
        if self.engine == 'psi4':
            self.input += f"#! Computation at {self.keywords[self.engine]['theory']}/{self.keywords[self.engine]['basis']} for {self.keywords[self.engine]['name']} \n\n"\
            f"memory {self.keywords[self.engine]['memory']}\n\n"\
            f"molecule {self.keywords[self.engine]['name']} {'{'}\n"\
            f"{self.keywords[self.engine]['charge']} {self.keywords[self.engine]['multiplicity']}\n"\
            f"{self.keywords[self.engine]['coords'].to_string(header=False, index=False)}\n"
            if self.keywords[self.engine]['PCM']: self.input += "symmetry c1\n"
            self.input += "}\n\n"
            if self.keywords[self.engine]['PCM']:
                self.input += f"set {'{'}\n  basis {self.keywords[self.engine]['basis']}\n"\
                "   scf_type pk\n"\
                "   pcm true\n"\
                "   pcm_scf_type total\n"\
                "   }\n\n"\
                "pcm = {\n"\
                "   Units = Angstrom\n"\
                "   Medium {\n"\
                "       SolverType = IEFPCM\n"\
                "       Solvent = Water\n"\
                "       }\n\n"\
                "   Cavity {\n"\
                "       RadiiSet = UFF\n"\
                "       Type = GePol\n"\
                "       Scaling = False\n"\
                "       Area = 0.3\n"\
                "       Mode = Implicit\n"\
                "       }\n"\
                "   }\n"
            else:
                self.input += f"set basis {self.keywords[self.engine]['basis']}\n"
            self.input += f"{self.keywords[self.engine]['calc_type']}('{self.keywords[self.engine]['theory']}')\n"
            if self.keywords[self.engine]['freq']:
               self.input += f"{self.keywords[self.engine]['theory']}_e, {self.keywords[self.engine]['theory']}_wfn = frequencies('{self.keywords[self.engine]['theory']}', return_wfn=True, dertype=1)\n"
            self.input += f"with open('{self.keywords[self.engine]['name']}_psi4out.xyz', 'w') as f:\n"\
            f"\tf.write('{len(self.keywords[self.engine]['coords'])} ' )\n"\
            f"\tf.write('%.12f\\n' % {self.keywords[self.engine]['theory']}_e)\n"\
            f"\tf.write({self.keywords[self.engine]['name']}.save_string_xyz())"

            self.jobsh += "#!/bin/bash\n"\
            f"#SBATCH --partition {self.keywords[self.engine]['partition']}\n"\
            f"#SBATCH --output myjob_{self.keywords[self.engine]['name']}.out\n"\
            f"#SBATCH --error myjob_{self.keywords[self.engine]['name']}.err\n"\
            f"#SBATCH --ntasks {self.keywords[self.engine]['ntasks']}\n"\
            f"#SBATCH --cpus-per-task {self.keywords[self.engine]['cpus-per-task']}\n"\
            f"#SBATCH --time {self.keywords[self.engine]['time']}\n"\
            f"#SBATCH --job-name={self.keywords[self.engine]['name']}\n"\
            f"#SBATCH --nodes {self.keywords[self.engine]['nodes']}\n"\
            f"#SBATCH --nice={self.keywords[self.engine]['nice']}\n"\
            f"#SBATCH --gpus={self.keywords[self.engine]['gpus']}\n"
            if self.keywords[self.engine]['mail_user']:
                self.jobsh += f"#SBATCH --mail-user={self.keywords[self.engine]['mail_user']}\n"
            if self.keywords[self.engine]['exclude']:
                self.jobsh +=f"#SBATCH --exclude={self.keywords[self.engine]['exclude']}\n"

            if self.machine == 'gwdg':
                self.jobsh += f"#SBATCH -A all\n"\
                "#SBATCH -C scratch\n"            
            
            self.jobsh +="\n# This block is for the execution of the program\n"
            if self.machine == 'smaug':
                self.jobsh += "source /home/users/all-jh/opt/miniconda3/etc/profile.d/conda.sh #It is needed to use the conda activate command\n"\
                "conda activate htmd\n\n"\
                "# Creating local scratch folder for the user on the computing node.\n"\
                "MY_TEMP_DIR=\"$(mktemp -d /localdisk/psi4_${SLURM_JOBID}_$(date +%Y-%m-%d-%H-%M-%S)-XXXXXXXXXX)\"\n"
            elif self.machine == 'gwdg':
                self.jobsh += "module load anaconda3/2020.07\n"\
                "source /opt/sw/rev/20.12/haswell/gcc-9.3.0/anaconda3-2020.07-slbv7z/etc/profile.d/conda.sh #It is needed to use the conda activate command\n"\
                "conda activate htmd\n\n"\
                "# Creating local scratch folder for the user on the computing node.\n"\
                "MY_TEMP_DIR=\"$(mktemp -d /scratch/users/${USER}/psi4_${SLURM_JOBID}_$(date +%Y-%m-%d-%H-%M-%S)-XXXXXXXXXX)\"\n"
            else:
                pass

            self.jobsh += "# make sure the temp dir was created\n"\
            "[ -d $MY_TEMP_DIR ] || exit 1\n"\
            "# make sure to delete the directory on (erroneous) exit\n"\
            "trap 'rm -rf -- \"$MY_TEMP_DIR\"' EXIT\n"\
            "# export PSI_SCRATCH env variable\n"\
            "export PSI_SCRATCH=$MY_TEMP_DIR\n\n"\
            "# This block is echoing some SLURM variables\n"\
            "echo \"Job execution start: $(date)\"\n"\
            "echo \"Jobid = $SLURM_JOBID\"\n"\
            "echo \"Host = $SLURM_JOB_NODELIST\"\n"\
            "echo \"Jobname = $SLURM_JOB_NAME\"\n"\
            "echo \"Subcwd = $SLURM_SUBMIT_DIR\"\n"\
            "echo \"SLURM_CPUS_PER_TASK = $SLURM_CPUS_PER_TASK\"\n"\
            "echo \"SLURM_CPUS_ON_NODE = $SLURM_CPUS_ON_NODE\"\n"\
            "echo \"The temp file used was: $MY_TEMP_DIR\"\n\n"\
            f"psi4 -i ${{SLURM_JOB_NAME}}.in -o ${{SLURM_JOB_NAME}}.out -n {self.keywords[self.engine]['ntasks'] * self.keywords[self.engine]['cpus-per-task']} #run command, in, out and number of threads to be used\n\n"\
            "#Deliting the scratch file. Normal delete on normal exit\n"\
            "rm -rf $MY_TEMP_DIR\n\n"\
            "echo \"Job execution end: $(date)\""
 
        elif self.engine == 'orca':
            #!! This is not implemented for GWDG
            self.input += f"#! Computation at {self.keywords[self.engine]['theory']}/{self.keywords[self.engine]['basis']} for {self.keywords[self.engine]['name']}\n"\
            f"%pal nprocs {self.keywords[self.engine]['ntasks']} end\n"
            if self.keywords[self.engine]['freq_type']:
                self.input += f"%maxcore {self.keywords[self.engine]['maxcore']}\n\n"\
                f"! {self.keywords[self.engine]['SCF_details']}\n"\
                f"! {self.keywords[self.engine]['freq_type']}\n"
            self.input += f"! {self.keywords[self.engine]['theory']} {self.keywords[self.engine]['basis']}\n"\
            f"! {self.keywords[self.engine]['calc_type']}\n"
            if self.keywords[self.engine]['SMD']:
                self.input += "%cpcm\n"\
                "   smd true\n"\
                "   SMDsolvent \"water\"\n"\
                "end\n"
            self.input += f"\n* xyz {self.keywords[self.engine]['charge']} {self.keywords[self.engine]['multiplicity']}\n"\
            f"! {self.keywords[self.engine]['coords'].to_string(header=False, index=False)}\n*"

            self.jobsh += "#!/bin/bash\n"\
            f"#SBATCH --partition {self.keywords[self.engine]['partition']}\n"\
            f"#SBATCH --output myjob_{self.keywords[self.engine]['name']}.out\n"\
            f"#SBATCH --error myjob_{self.keywords[self.engine]['name']}.err\n"\
            f"#SBATCH --ntasks {self.keywords[self.engine]['ntasks']}\n"\
            f"#SBATCH --cpus-per-task {self.keywords[self.engine]['cpus-per-task']}\n"\
            f"#SBATCH --time {self.keywords[self.engine]['time']}\n"\
            f"#SBATCH --job-name={self.keywords[self.engine]['name']}\n"\
            f"#SBATCH --nodes {self.keywords[self.engine]['nodes']}\n"\
            f"#SBATCH --nice={self.keywords[self.engine]['nice']}\n"\
            f"#SBATCH --gpus={self.keywords[self.engine]['gpus']}\n"
            if self.keywords[self.engine]['mail_user']:
                self.jobsh += f"#SBATCH --mail-user={self.keywords[self.engine]['mail_user']}\n"
            if self.keywords[self.engine]['exclude']:
                self.jobsh +=f"#SBATCH --exclude={self.keywords[self.engine]['exclude']}\n"
            
            self.jobsh +="\n# This block is for the execution of the program\n"\
            "# Setting OPENMPI paths here:\n"\
            "export PATH=\"/data/shared/opt/ORCA/openmpi314/bin:$PATH\"\n"\
            "export LD_LIBRARY_PATH=\"$LD_LIBRARY_PATH:/data/shared/opt/ORCA/openmpi314/lib\"\n\n"\
            "# Here giving the path to the ORCA binaries and giving communication protocol\n"\
            "export PATH=\"/data/shared/opt/ORCA/orca_4_2_1_linux_x86-64_openmpi314:$PATH\"\n"\
            "export LD_LIBRARY_PATH=\"/data/shared/opt/ORCA/orca_4_2_1_linux_x86-64_openmpi314:$LD_LIBRARY_PATH\"\n\n"\
            "# Creating local scratch folder for the user on the computing node.\n"\
            "MY_TEMP_DIR=\"$(mktemp -d /localdisk/psi4_${SLURM_JOBID}_$(date +%Y-%m-%d-%H-%M-%S)-XXXXXXXXXX)\"\n"\
            "# make sure the temp dir was created\n"\
            "[ -d $MY_TEMP_DIR ] || exit 1\n"\
            "# make sure to delete the directory on (erroneous) exit\n"\
            "trap 'rm -rf -- \"$MY_TEMP_DIR\"' EXIT\n\n"\
            "# Copy only the necessary stuff in submit directory to scratch directory. Add more here if needed.\n"\
            "cp $SLURM_SUBMIT_DIR/*.inp $MY_TEMP_DIR/\n"\
            "#cp $SLURM_SUBMIT_DIR/*.gbw $MY_TEMP_DIR/\n"\
            "#cp $SLURM_SUBMIT_DIR/*.hess $MY_TEMP_DIR/\n\n"\
            "# Creating nodefile in scratch\n"\
            "echo $SLURM_NODELIST > $MY_TEMP_DIR/${SLURM_JOB_NAME}.nodes\n\n"\
            "# This block is echoing some SLURM variables\n"\
            "echo \"Job execution start: $(date)\"\n"\
            "echo \"Jobid = $SLURM_JOBID\"\n"\
            "echo \"Host = $SLURM_JOB_NODELIST\"\n"\
            "echo \"Jobname = $SLURM_JOB_NAME\"\n"\
            "echo \"Subcwd = $SLURM_SUBMIT_DIR\"\n"\
            "echo \"SLURM_CPUS_PER_TASK = $SLURM_CPUS_PER_TASK\"\n"\
            "echo \"SLURM_CPUS_ON_NODE = $SLURM_CPUS_ON_NODE\"\n"\
            "echo \"The temp file used was: $MY_TEMP_DIR\"\n\n"\
            f"$(which orca) ${{SLURM_JOB_NAME}}.inp > ${{SLURM_JOB_NAME}}.log 2>&1 # run command\n\n"\
            "# ORCA has finished here. Now copy important stuff back (xyz files, GBW files etc.). Add more here if needed.\n"\
            "cp $MY_TEMP_DIR/*.gbw $SLURM_SUBMIT_DIR\n"\
            "cp $MY_TEMP_DIR/*.xyz $SLURM_SUBMIT_DIR\n"\
            "cp $MY_TEMP_DIR/*.opt $SLURM_SUBMIT_DIR\n"\
            "cp $MY_TEMP_DIR/*.hess $SLURM_SUBMIT_DIR\n"\
            "cp $MY_TEMP_DIR/*_property.txt $SLURM_SUBMIT_DIR\n"\
            "cp $MY_TEMP_DIR/*.log $SLURM_SUBMIT_DIR\n\n"\
            "#Deleting the scratch file. Normal delete on normal exit\n"\
            "rm -rf $MY_TEMP_DIR\n\n"\
            "echo \"Job execution end: $(date)\""
        elif self.engine == 'gaussian':
            # Right now we don't have access to gaussian, so, it will not be coded.
            # But is more than the same, configurate the .gjf file and the corresponded .sh
            raise Exception(f'The engine {self.engine} is not coded!')
        else:
            raise Exception(f'The engine {self.engine} is not coded!')
    def write(self, outpath, attribute):
        with open(outpath, 'w') as f:
            if attribute == 'input':
                f.write(self.input)
            elif attribute == 'jobsh':
                f.write(self.jobsh)
            else:
                raise ValueError(f"{attribute} is not a correct attribute. Must be: 'input' or 'jobsh'")

class PARAM:
    def __init__(self, machine = 'smaug', **keywords):
        self.machine = machine
        self.default_partition = {'smaug': 'deflt', 'gwdg': 'medium'}
        self.default_exclude = {'smaug': "fang[1,11-50]", 'gwdg': None}
        self.default_ntasks = {'smaug': 12, 'gwdg': 24}
        self.default_keywords = {
                #Parameters for input
                'theory': 'MP2',
                'basis': 'aug-cc-pVTZ',
                'memory': 14500,
                'name': 'mol',
                'environment': 'vacuum',
                'charge': 0,
                'coords': pd.DataFrame(),
                'min-type': None,
                'charge-type': 'ESP',
                'scan-type': 'mm',
                'dihed-fit-type': 'iterative',
                'dihed-num-iterations': 100,
                'outdir': 'outdir',
                #Parameters for jobsh
                'time': '2-00:00',
                'ntasks': self.default_ntasks[self.machine],
                'cpus-per-task': 1,
                'partition': self.default_partition[self.machine],
                'nodes': 1,
                'nice': 0,
                'gpus': 0,
                'mail_user': None,
                'exclude': self.default_partition[self.machine]}
        self.keywords = copy.deepcopy(self.default_keywords)
        for key in keywords:
            self.keywords[key] = keywords[key]
        self.jobsh = ''
        self.parse()
   
    def parse(self):
        self.jobsh += "#!/bin/bash\n"\
        f"#SBATCH --partition {self.keywords['partition']}\n"\
        f"#SBATCH --output myjob_{self.keywords['name']}.out\n"\
        f"#SBATCH --error myjob_{self.keywords['name']}.err\n"\
        f"#SBATCH --ntasks {self.keywords['ntasks']}\n"\
        f"#SBATCH --cpus-per-task {self.keywords['cpus-per-task']}\n"\
        f"#SBATCH --time {self.keywords['time']}\n"\
        f"#SBATCH --job-name={self.keywords['name']}\n"\
        f"#SBATCH --nodes {self.keywords['nodes']}\n"\
        f"#SBATCH --nice={self.keywords['nice']}\n"\
        f"#SBATCH --gpus={self.keywords['gpus']}\n"
        if self.keywords['mail_user']:
            self.jobsh += f"#SBATCH --mail-user={self.keywords['mail_user']}\n"
        if self.keywords['exclude']:
            self.jobsh +=f"#SBATCH --exclude={self.keywords['exclude']}\n"
    
        if self.machine == 'gwdg':
            self.jobsh += f"#SBATCH -A all\n"\
            "#SBATCH -C scratch\n"            
        
        self.jobsh +="\n# This block is for the execution of the program\n"
        if self.machine == 'smaug':
            self.jobsh += "source /home/users/all-jh/opt/miniconda3/etc/profile.d/conda.sh #It is needed to use the conda activate command\n"\
            "conda activate htmd\n\n"\
            "# Creating local scratch folder for the user on the computing node.\n"\
            "MY_TEMP_DIR=\"$(mktemp -d /localdisk/psi4_${SLURM_JOBID}_$(date +%Y-%m-%d-%H-%M-%S)-XXXXXXXXXX)\"\n"
        elif self.machine == 'gwdg':
            self.jobsh += "module load anaconda3/2020.07\n"\
            "source /opt/sw/rev/20.12/haswell/gcc-9.3.0/anaconda3-2020.07-slbv7z/etc/profile.d/conda.sh #It is needed to use the conda activate command\n"\
            "conda activate htmd\n\n"\
            "# Creating local scratch folder for the user on the computing node.\n"\
            "MY_TEMP_DIR=\"$(mktemp -d /scratch/users/${USER}/psi4_${SLURM_JOBID}_$(date +%Y-%m-%d-%H-%M-%S)-XXXXXXXXXX)\"\n"
        else:
            pass

        self.jobsh += "# make sure the temp dir was created\n"\
        "[ -d $MY_TEMP_DIR ] || exit 1\n"\
        "# make sure to delete the directory on (erroneous) exit\n"\
        "trap 'rm -rf -- \"$MY_TEMP_DIR\"' EXIT\n"\
        "# export PSI_SCRATCH env variable\n"\
        "export PSI_SCRATCH=$MY_TEMP_DIR\n\n"\
        "# This block is echoing some SLURM variables\n"\
        "echo \"Job execution start: $(date)\"\n"\
        "echo \"Jobid = $SLURM_JOBID\"\n"\
        "echo \"Host = $SLURM_JOB_NODELIST\"\n"\
        "echo \"Jobname = $SLURM_JOB_NAME\"\n"\
        "echo \"Subcwd = $SLURM_SUBMIT_DIR\"\n"\
        "echo \"SLURM_CPUS_PER_TASK = $SLURM_CPUS_PER_TASK\"\n"\
        "echo \"SLURM_CPUS_ON_NODE = $SLURM_CPUS_ON_NODE\"\n"\
        "echo \"The temp file used was: $MY_TEMP_DIR\"\n\n"\
        "#run paremeterize\n"\
        f"parameterize ${{SLURM_JOB_NAME}}.mol2 --memory {self.keywords['memory']} --ncpus {self.keywords['ntasks'] * self.keywords['cpus-per-task']} --theory {self.keywords['theory']} --basis {self.keywords['basis']} --environment {self.keywords['environment']} --charge {self.keywords['charge']} --min-type {self.keywords['min-type']} --charge-type {self.keywords['charge-type']} --scan-type {self.keywords['scan-type']} --dihed-fit-type {self.keywords['dihed-fit-type']} --dihed-num-iterations {self.keywords['dihed-num-iterations']} --outdir {self.keywords['outdir']} "
        "#Deleting the scratch file. Normal delete on normal exit\n"\
        "rm -rf $MY_TEMP_DIR\n\n"\
        "echo \"Job execution end: $(date)\""

    def write(self, outpath = 'job.sh'):
        with open(outpath, 'w') as f:
            f.write(self.jobsh)



class CONTINUE:
    def __init__(self, elapsed_paths, engine = 'psi4', machine = 'smaug', **keywords):
        self.elapsed_paths = elapsed_paths
        self.engine = engine
        self.machine = machine
        self.default_partition = {'smaug': 'deflt', 'gwdg': 'medium'}
        self.default_exclude = {'smaug': "fang[1,11-50]", 'gwdg': None}
        self.default_ntasks = {'smaug': 12, 'gwdg': 24}
        self.default_keywords = {
            'psi4':{
                #Parameters for jobsh
                'time': '2-00:00',
                'ntasks': self.default_ntasks[self.machine],
                'cpus-per-task': 1,
                'partition': self.default_partition[self.machine],
                'nodes': 1,
                'nice': 0,
                'gpus': 0,
                'mail_user': None,
                'exclude': self.default_exclude[self.machine],
                },
            'orca':{
                #Parameters for jobsh
                'time': '2-00:00',
                'ntasks': 6,
                'cpus-per-task': 2,
                'partition': 'deflt',
                'nodes': 1,
                'nice': 0,
                'gpus': 0,
                'mail_user': None,
                'exclude': "fang[1,11-50]" #None
                },
            'gaussian':{
                'theory': 7,
                'basis': 2,
                'memory':4,
                'name':5,
                'charge':0,
                'multiplicity':1,
                'PCM': False,
                'coords': None,
                'calc_type':5,
                'freq': False}
        }
        self.keywords = copy.deepcopy(self.default_keywords)
        for key in keywords:
            self.keywords[engine][key] = keywords[key]
        self.jobsh = ''
        self.parse()
    
    def parse(self):
        if self.engine == 'psi4':

            self.jobsh += "#!/bin/bash\n"\
            f"#SBATCH --partition {self.keywords[self.engine]['partition']}\n"\
            f"#SBATCH --output myjob_{self.keywords[self.engine]['name']}.out\n"\
            f"#SBATCH --error myjob_{self.keywords[self.engine]['name']}.err\n"\
            f"#SBATCH --ntasks {self.keywords[self.engine]['ntasks']}\n"\
            f"#SBATCH --cpus-per-task {self.keywords[self.engine]['cpus-per-task']}\n"\
            f"#SBATCH --time {self.keywords[self.engine]['time']}\n"\
            f"#SBATCH --job-name={self.keywords[self.engine]['name']}\n"\
            f"#SBATCH --nodes {self.keywords[self.engine]['nodes']}\n"\
            f"#SBATCH --nice={self.keywords[self.engine]['nice']}\n"\
            f"#SBATCH --gpus={self.keywords[self.engine]['gpus']}\n"
            if self.keywords[self.engine]['mail_user']:
                self.jobsh += f"#SBATCH --mail-user={self.keywords[self.engine]['mail_user']}\n"
            if self.keywords[self.engine]['exclude']:
                self.jobsh +=f"#SBATCH --exclude={self.keywords[self.engine]['exclude']}\n"

            if self.machine == 'gwdg':
                self.jobsh += f"#SBATCH -A all\n"\
                "#SBATCH -C scratch\n"            
            
            self.jobsh +="\n# This block is for the execution of the program\n"
            if self.machine == 'smaug':
                self.jobsh += "source /home/users/all-jh/opt/miniconda3/etc/profile.d/conda.sh #It is needed to use the conda activate command\n"\
                "conda activate htmd\n\n"\
                "# Creating local scratch folder for the user on the computing node.\n"\
                "MY_TEMP_DIR=\"$(mktemp -d /localdisk/psi4_${SLURM_JOBID}_$(date +%Y-%m-%d-%H-%M-%S)-XXXXXXXXXX)\"\n"
            elif self.machine == 'gwdg':
                self.jobsh += "module load anaconda3/2020.07\n"\
                "source /opt/sw/rev/20.12/haswell/gcc-9.3.0/anaconda3-2020.07-slbv7z/etc/profile.d/conda.sh #It is needed to use the conda activate command\n"\
                "conda activate htmd\n\n"\
                "# Creating local scratch folder for the user on the computing node.\n"\
                "MY_TEMP_DIR=\"$(mktemp -d /scratch/users/${USER}/psi4_${SLURM_JOBID}_$(date +%Y-%m-%d-%H-%M-%S)-XXXXXXXXXX)\"\n"
            else:
                pass

            self.jobsh += "# make sure the temp dir was created\n"\
            "[ -d $MY_TEMP_DIR ] || exit 1\n"\
            "# make sure to delete the directory on (erroneous) exit\n"\
            "trap 'rm -rf -- \"$MY_TEMP_DIR\"' EXIT\n"\
            "# export PSI_SCRATCH env variable\n"\
            "export PSI_SCRATCH=$MY_TEMP_DIR\n\n"\
            "# This block is echoing some SLURM variables\n"\
            "echo \"Job execution start: $(date)\"\n"\
            "echo \"Jobid = $SLURM_JOBID\"\n"\
            "echo \"Host = $SLURM_JOB_NODELIST\"\n"\
            "echo \"Jobname = $SLURM_JOB_NAME\"\n"\
            "echo \"Subcwd = $SLURM_SUBMIT_DIR\"\n"\
            "echo \"SLURM_CPUS_PER_TASK = $SLURM_CPUS_PER_TASK\"\n"\
            "echo \"SLURM_CPUS_ON_NODE = $SLURM_CPUS_ON_NODE\"\n"\
            "echo \"The temp file used was: $MY_TEMP_DIR\"\n\n"
            for path_and_file in self.elapsed_paths:
                path, file = os.path.split(path_and_file)
                name = file.split('.')[0]
                self.jobsh += f"psi4 -i {os.path.join(path, name)}.in -o {os.path.join(path, name)}.out &> {os.path.join(path, name)}.log\n"\

            f"psi4 -i ${{SLURM_JOB_NAME}}.in -o ${{SLURM_JOB_NAME}}.out -n {self.keywords[self.engine]['ntasks'] * self.keywords[self.engine]['cpus-per-task']} #run command, in, out and number of threads to be used\n\n"\
            "#Deliting the scratch file. Normal delete on normal exit\n"\
            "rm -rf $MY_TEMP_DIR\n\n"\
            "echo \"Job execution end: $(date)\""
 
        elif self.engine == 'orca':
            raise Exception(f'The engine {self.engine} is not coded!')
            # Not implemented
            self.jobsh += "#!/bin/bash\n"\
            f"#SBATCH --partition {self.keywords[self.engine]['partition']}\n"\
            f"#SBATCH --output myjob_{self.keywords[self.engine]['name']}.out\n"\
            f"#SBATCH --error myjob_{self.keywords[self.engine]['name']}.err\n"\
            f"#SBATCH --ntasks {self.keywords[self.engine]['ntasks']}\n"\
            f"#SBATCH --cpus-per-task {self.keywords[self.engine]['cpus-per-task']}\n"\
            f"#SBATCH --time {self.keywords[self.engine]['time']}\n"\
            f"#SBATCH --job-name={self.keywords[self.engine]['name']}\n"\
            f"#SBATCH --nodes {self.keywords[self.engine]['nodes']}\n"\
            f"#SBATCH --nice={self.keywords[self.engine]['nice']}\n"\
            f"#SBATCH --gpus={self.keywords[self.engine]['gpus']}\n"
            if self.keywords[self.engine]['mail_user']:
                self.jobsh += f"#SBATCH --mail-user={self.keywords[self.engine]['mail_user']}\n"
            if self.keywords[self.engine]['exclude']:
                self.jobsh +=f"#SBATCH --exclude={self.keywords[self.engine]['exclude']}\n"
            
            self.jobsh +="\n# This block is for the execution of the program\n"\
            "# Setting OPENMPI paths here:\n"\
            "export PATH=\"/data/shared/opt/ORCA/openmpi314/bin:$PATH\"\n"\
            "export LD_LIBRARY_PATH=\"$LD_LIBRARY_PATH:/data/shared/opt/ORCA/openmpi314/lib\"\n\n"\
            "# Here giving the path to the ORCA binaries and giving communication protocol\n"\
            "export PATH=\"/data/shared/opt/ORCA/orca_4_2_1_linux_x86-64_openmpi314:$PATH\"\n"\
            "export LD_LIBRARY_PATH=\"/data/shared/opt/ORCA/orca_4_2_1_linux_x86-64_openmpi314:$LD_LIBRARY_PATH\"\n\n"\
            "# Creating local scratch folder for the user on the computing node.\n"\
            "MY_TEMP_DIR=\"$(mktemp -d /localdisk/psi4_${SLURM_JOBID}_$(date +%Y-%m-%d-%H-%M-%S)-XXXXXXXXXX)\"\n"\
            "# make sure the temp dir was created\n"\
            "[ -d $MY_TEMP_DIR ] || exit 1\n"\
            "# make sure to delete the directory on (erroneous) exit\n"\
            "trap 'rm -rf -- \"$MY_TEMP_DIR\"' EXIT\n\n"\
            "# Copy only the necessary stuff in submit directory to scratch directory. Add more here if needed.\n"\
            "cp $SLURM_SUBMIT_DIR/*.inp $MY_TEMP_DIR/\n"\
            "#cp $SLURM_SUBMIT_DIR/*.gbw $MY_TEMP_DIR/\n"\
            "#cp $SLURM_SUBMIT_DIR/*.hess $MY_TEMP_DIR/\n\n"\
            "# Creating nodefile in scratch\n"\
            "echo $SLURM_NODELIST > $MY_TEMP_DIR/${SLURM_JOB_NAME}.nodes\n\n"\
            "# This block is echoing some SLURM variables\n"\
            "echo \"Job execution start: $(date)\"\n"\
            "echo \"Jobid = $SLURM_JOBID\"\n"\
            "echo \"Host = $SLURM_JOB_NODELIST\"\n"\
            "echo \"Jobname = $SLURM_JOB_NAME\"\n"\
            "echo \"Subcwd = $SLURM_SUBMIT_DIR\"\n"\
            "echo \"SLURM_CPUS_PER_TASK = $SLURM_CPUS_PER_TASK\"\n"\
            "echo \"SLURM_CPUS_ON_NODE = $SLURM_CPUS_ON_NODE\"\n"\
            "echo \"The temp file used was: $MY_TEMP_DIR\"\n\n"\
            f"$(which orca) ${{SLURM_JOB_NAME}}.inp > ${{SLURM_JOB_NAME}}.log 2>&1 # run command\n\n"\
            "# ORCA has finished here. Now copy important stuff back (xyz files, GBW files etc.). Add more here if needed.\n"\
            "cp $MY_TEMP_DIR/*.gbw $SLURM_SUBMIT_DIR\n"\
            "cp $MY_TEMP_DIR/*.xyz $SLURM_SUBMIT_DIR\n"\
            "cp $MY_TEMP_DIR/*.opt $SLURM_SUBMIT_DIR\n"\
            "cp $MY_TEMP_DIR/*.hess $SLURM_SUBMIT_DIR\n"\
            "cp $MY_TEMP_DIR/*_property.txt $SLURM_SUBMIT_DIR\n"\
            "cp $MY_TEMP_DIR/*.log $SLURM_SUBMIT_DIR\n\n"\
            "#Deleting the scratch file. Normal delete on normal exit\n"\
            "rm -rf $MY_TEMP_DIR\n\n"\
            "echo \"Job execution end: $(date)\""
        elif self.engine == 'gaussian':
            # Right now we don't have access to gaussian, so, it will not be coded.
            # But is more than the same, configurate the .gjf file and the corresponded .sh
            raise Exception(f'The engine {self.engine} is not coded!')
        else:
            raise Exception(f'The engine {self.engine} is not coded!')
    def write(self, outpath, attribute):
        with open(outpath, 'w') as f:
            if attribute == 'input':
                f.write(self.input)
            elif attribute == 'jobsh':
                f.write(self.jobsh)
            else:
                raise ValueError(f"{attribute} is not a correct attribute. Must be: 'input' or 'jobsh'")



if __name__ == '__main__':

    s = CONTINUE(['./psi4.in', 'klusef/fsdg/oi/psi4.in'],machine = 'gwdg', name = 'elapse1')
    #s = INPUT('psi4',machine = 'gwdg', name =  'imi')
    print(s.jobsh)