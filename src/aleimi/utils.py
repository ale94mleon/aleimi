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
import os, subprocess, time, inspect

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

        

def run(command:str, shell:bool = True, executable:str = '/bin/bash', Popen:bool = False):
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
def mopac(mop:str):
    """A simple wrapper around MOPAC

    Parameters
    ----------
    mop : str
        The path where the .mop file is located
    """
    print(f"Mopac is running ...")
    run(f"mopac {mop}  > /dev/null 2>&1")
    print("Done!")

def makedirs(path):
    os.makedirs(path,exist_ok=True)

def KbT(absolute_temperature):
    """Return the value of Kb*T in kJ/mol

    Args:
        absolute_temperature (float): The absolute temperature in kelvin
    """
    Kb = 8.314462618E-3 #kJ/(mol⋅K) (kNA)
    return absolute_temperature*Kb
def get_default_kwargs(func):
    signature = inspect.signature(func)
    return {
        k: v.default
        for k, v in signature.parameters.items()
        if v.default is not inspect.Parameter.empty
    }

if __name__ == '__main__':...

  