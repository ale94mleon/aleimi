#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
Created on    : Fri Jun 25 21:00:10 2021
Author        : Alejandro Martínez León
Mail          : [alejandro.martinezleon@uni-saarland.de, ale94mleon@gmail.com]
Affiliation   : Jochen Hub's Biophysics Group
Affiliation   : Faculty of NS, University of Saarland, Saarbrücken, Germany
===============================================================================
DESCRIPTION   :
DEPENDENCIES  :
===============================================================================
"""
import os, subprocess, shutil, datetime, tempfile, time, tqdm, inspect
from matplotlib.pyplot import get_figlabels
import numpy as np
import multiprocessing as mp
from glob import glob
from mdynamic.tools import xvg
#=======================================================================================

#                          Miscellanea tools

#=======================================================================================

def zerolistmaker(n):
    return [0]*n
def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
def aovec(vectors, round = None):
    """Return the average oriented vector
    Calculate all the angles respect to the cartessina axis and average the angles,andthen return the unitary vector
    ang(vec, x) = mean(ang(vectors_i, x))
    ang(vec, y) = mean(ang(vectors_i, y))
    ang(vec, z) = mean(ang(vectors_i, z))
    vec = (cos(ang(vec, x); cos(ang(vec, y); cos(ang(vec, z))
    Args:
        vectors ([type]): [description]
        round ([type]): [description]
    Returns:
        [type]: [description]
    """
    vectors = [np.asarray(v) for v in vectors]
    #Dimension of the vectors, how many components
    vec = np.empty(vectors[0].shape)
    for i in range(vectors[0].shape[0]):
        zero_vector = np.zeros(vectors[0].shape)
        zero_vector[i] = 1
        vec[i] = np.cos(np.mean([angle_between(vector, zero_vector) for vector in vectors]))

    vec = unit_vector(vec)
    if round: vec = np.round(vec, round)
    return vec



def timeit(method):
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        if 'log_time' in kw:
            name = kw.get('log_name', method.__name__.upper())
            kw['log_time'][name] = int((te - ts) * 1000)
        else:
            print('%r  %2.2f ms' % \
                  (method.__name__, (te - ts) * 1000))
        return result
    return timed

#=======================================================================================

#                          Tools for execution

#=======================================================================================
def multi_run(commands, nPar, shell = True, executable = '/bin/bash'):
    """This will run as many runs as nPar.

    Args:
        commands (list): A list of string to be run in the specified shell.
        nPar (int): How many processes are running simultaneously
        shell (bool, optional): belongs to run(). Defaults to True.
        executable (str, optional): belongs to tun(). Defaults to '/bin/bash'.
        Popen (bool, optional): belongs to run(). Defaults to False.
    """

    tmp_cmds = []
    for cmd in commands:
        tmp_cmds.append(run(cmd, shell = shell, executable = executable, Popen = True))
        if len(tmp_cmds) >= nPar:
            print(f"{len(tmp_cmds)} running, now waiting...")
            for tmp_cmd in tmp_cmds:
                tmp_cmd.wait()
            tmp_cmds = []
        

def run(command, shell = True, executable = '/bin/bash', Popen = False):
    #Here I could make some modification in order that detect the operator system
    #NAd make the command compatible with the opertor system
    #the fucntion eval could be an option if some modifcation to the variable command 
    #need to be done.... SOme fligth ideas...

    if Popen:
        #In this case you could acces the pid as: run.pid
        process = subprocess.Popen(command, shell = shell, executable = executable)
    else:
        process = subprocess.run(command, shell = shell, executable = executable)
    return process
@timeit
def mopac(mop, mopacExecutablePath = '/opt/mopac/MOPAC2016.exe'):
    print(f"Mopac is running ...")
    run(f"echo | {mopacExecutablePath} {mop}  > /dev/null 2>&1")
    print("Done!")


def checkrun():
    """Check for the running process on the cluster.

    Returns:
        list of integers: The integers that identify the running process on the cluster.
    """
    tmp = tempfile.NamedTemporaryFile()
    run(f"squeue -u $USER --format=%.i > {tmp.name}")
    with open(tmp.name, "r") as f:
        lines = f.readlines()
    return [int(item.strip()) for  item in lines[1:]]#The first line is the string "JOBID"

def job_launch(shell="sbatch", script_name = "job.sh"):
    """
    

    Parameters
    ----------
    shell : TYPE, optional
        DESCRIPTION. The default is "sbatch".
    script_name : TYPE, optional
        DESCRIPTION. The default is "job.sh".
        If a regular expresion is provided
        then the function will execute the first ordered alphabetically. E.g:
            job.* was provided and there job.sh and job.bash. Then it will use job.bash.


    Returns
    -------
    JOBIDs : list of integers.
        DESCRIPTION. The ID of the launch in case of sbatch was used as shell

    """
    cwd = os.getcwd()
    JOBIDs = []
    for (root, dirs, files) in list(os.walk(cwd)):
        os.chdir(root)
        if glob(script_name):
            script_name = sorted(glob(script_name))[0]
            JOBID_tmp_file = tempfile.NamedTemporaryFile()
            if shell == "sbatch":
                run(f"{shell} --parsable {script_name}>{JOBID_tmp_file.name}")
                with open(JOBID_tmp_file.name, "r") as f:
                    JOBIDs.append(int(f.readline()))
            else:
                run(f"{shell} {script_name}>{JOBID_tmp_file.name}")
                try:
                    with open(JOBID_tmp_file.name, "r") as f:
                        JOBIDs.append(int(f.readline()))
                except:
                    pass

    os.chdir(cwd)
    return JOBIDs


def backoff(file_path):
    """
    

    Parameters
    ----------
    file : TYPE string
        DESCRIPTION: The name or the path f     or the specific file

    Returns
    -------
    None.
    If the file already exist. it will made a back up to ./#{file}.{str(i)}#,
    Where i is an integer.
    """
    basname = os.path.basename(file_path)
    dirname = os.path.dirname(file_path)
    if os.path.exists(file_path):
        new_basname = basname
        i = 1
        while(os.path.exists(os.path.join(dirname, new_basname))):
            new_basname = f"./#{basname}.{str(i)}#"
            i += 1
        print(f"Back Off! I just backed up {file_path} to {os.path.join(dirname, new_basname)}")
        shutil.copy2(file_path, os.path.join(dirname, new_basname))
#=======================================================================================

#                          Tools for working with files

#=======================================================================================
def makedirs(path):
    os.makedirs(path,exist_ok=True)

def rm(pattern, r = False):
    """
    

    Parameters
    ----------
    patterns : string
        input-like Unix rm
    r : TYPE, bool
        DESCRIPTION. The default is False.
        If True delete also directories
    Returns
    -------
    None.

    """
    if glob(pattern):
        for p in glob(pattern):

            if r:
                if os.path.isdir(p):
                    shutil.rmtree(p)
                elif os.path.isfile(p):
                    os.remove(p)       
            else:
                if os.path.isfile(p):
                    os.remove(p)
                elif os.path.isdir(p):
                    print(f"The directory '{p}' was not deleted, set the flag r to True")
    else:
        
        print(f"rm: '{pattern}' doesn't exist")

                
def cp(src, dest, r = False):
    """
    
    This function makes use of the possible multiple CPU of the machine.
    Parameters
    ----------
    src : TYPE: string
        DESCRIPTION:
            Source Path to a directory or a file. regular expresion are accepted. The librery glob is used for that.
    dest : TYPE: string
        DESCRIPTION:
            Destination Path
    r : TYPE, optional
        DESCRIPTION. The default is False:
            If True, Also directories will be copy. If not, and a directory was
            given as src, a Raise Exception will be printed
            Another Raise Exception will be printed if the destination path doesn't
            exist.

    Returns
    -------
    None.

    """
    src = os.path.abspath(src)
    dest = os.path.abspath(dest)
    #This scheme doesn't consider the possibilities of other item except files, dir or regualr expresions
    if glob(src) and os.path.exists(os.path.dirname(dest)):
    
        if os.path.isdir(src):
            if r == True:            
                
                basename_src = os.path.basename(src)
                for root, dirs, files in os.walk(src):
                    path2copy = root.split(src)[1]
                    try:
                        if path2copy[0] == "/":
                            path2copy = path2copy[1:]
                    except: pass
                    path_dest = os.path.join(dest,basename_src,path2copy)
                    
                    if dirs:
                        pool = mp.Pool(mp.cpu_count())
                        pool.map(makedirs, [os.path.join(path_dest, d) for d in dirs])
                        pool.close()
                    if files:
                        makedirs(path_dest)
                        pool = mp.Pool(mp.cpu_count())
                        
                        pool.starmap(shutil.copy2, [(os.path.join(root, f),path_dest) for f in files])
                        pool.close()
            else:
                print(f"If you need to copy directories, set the 'r' flag to True. {src} is a directory")
    
        elif os.path.isfile(src):
            shutil.copy2(src,dest)
            #if not os.path.exists(dest):
             #   print(f"The file '{src} was not copy to '{dest}' becasue doesn't exist or is not accesible")
        
        
        else:#Here we are in regular expresions
            to_copy = [file for file in glob(src) if os.path.isfile(file)]
            pool = mp.Pool(mp.cpu_count())
            pool.starmap(shutil.copy2, [(pattern,dest) for pattern in to_copy])
            pool.close()
        
    elif not glob(src) and not os.path.exists(os.path.dirname(dest)):
        print(f"cp: neither {src}, nor {dest} exist.")
    elif not glob(src):
        print(f"cp: The source file {src} doesn't exist")
    else:
        print(f"cp: cannot create regular file '{dest}': Not a directory")
            
        

def mv(src, dest, r = False):
    cp(src, dest, r = r)
    rm(src,r=r)

def list_if_dir(path = '.'):
    return [item for item in os.listdir(path) if os.path.isdir(os.path.join(path, item))]

def list_if_file(path = '.'):
    return [item for item in os.listdir(path) if os.path.isfile(os.path.join(path, item))]


def KbT(absolute_temperature):
    """Return the value of Kb*T in kJ/mol

    Args:
        absolute_temperature (float): The absolute temperature in kelvin
    """
    Kb = 8.314462618E-3 #kJ/(mol⋅K) (kNA)
    return absolute_temperature*Kb
def get_default_args(func):
    signature = inspect.signature(func)
    return {
        k: v.default
        for k, v in signature.parameters.items()
        if v.default is not inspect.Parameter.empty
    }

if __name__ == '__main__':  
    p = get_default_args(run)
    print(type(p['shell']))

  