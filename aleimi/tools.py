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
import os, subprocess, shutil, tempfile, time, inspect
import numpy as np
import multiprocessing as mp
from glob import glob




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
def mopac(mop):
    print(f"Mopac is running ...")
    run(f"mopac {mop}  > /dev/null 2>&1")
    print("Done!")


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

  