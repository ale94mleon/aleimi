#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, subprocess, time, inspect

def timeit(method):
    """Useful as decorator
    """
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

        

def run(command: str, shell: bool = True, executable: str = '/bin/bash'):
    """This function is just a useful wrapper around subprocess.run

    Parameters
    ----------
    command : str
        Any command to execute.
    shell : bool, optional
        keyword of ``subprocess.Popen`` and ``subprocess.Popen``, by default True
    executable : str, optional
        keyword of ``subprocess.Popen`` and ``subprocess.Popen``, by default '/bin/bash'

    Returns
    -------
    object
        The processes returned by Run.

    Raises
    ------
    RuntimeError
        In case of non-zero exit status on the provided command.
    """

    process = subprocess.run(command, shell=shell, executable=executable, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    returncode = process.returncode
    if returncode != 0:
        # print(f'Command {command} returned non-zero exit status {returncode}')
        raise RuntimeError(process.stderr)
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

def get_default_kwargs(func:object) -> dict:
    """Get the keywords arguments of a function

    Parameters
    ----------
    func : object
        The function to analysis

    Returns
    -------
    dict
        keyword:keyword_value
    """
    signature = inspect.signature(func)
    return {
        k: v.default
        for k, v in signature.parameters.items()
        if v.default is not inspect.Parameter.empty
    }

def ignoreLines(file_obj:object, num_of_lines_to_skip:int):
    """Ignore lines during reading of file

    Parameters
    ----------
    file_obj : object
        The file in open mode
    num_of_lines_to_skip : int
        Number of lines to skip
    """
    for _ in range(num_of_lines_to_skip): file_obj.readline()


if __name__ == '__main__':...

  