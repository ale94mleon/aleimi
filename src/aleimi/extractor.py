#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import os
from aleimi import templates, utils
from typing import List

# def extract(arc_file, boltzmann_file, energy_cut = 2, conformer_cut = None, mksh = True, mkdir = True):
# Energy cut is in kcal/mol


def extract(boltzmann_file: str, energy_cut: float = 2.0, conformer_cut: int = None) -> List[int]:
    """Extract the conformers based on the filters: energy_cut and/or conformer_cut

    Parameters
    ----------
    boltzmann_file : str
        The boltzmann file generated with `:meth: aleimi.boltzmann.main`
    energy_cut : float, optional
        Maximum difference in energy with respect the conformer with lowest energy, by default 2.0
    conformer_cut : int, optional
        Maximum number of conformers to export, by default None

    Returns
    -------
    List[int]
        The list of ``cell`` identifiers
    """

    df = pd.read_csv(boltzmann_file)
    indx_to_extract = []
    if energy_cut:
        df_subset = df[df.Emin_Ei >= -energy_cut]
        indx_to_extract += df_subset.cell.tolist()

    if conformer_cut:
        df_subset = df.iloc[:conformer_cut]
        indx_to_extract += df_subset.cell.tolist()

    if not energy_cut and not conformer_cut:
        df_subset = df.iloc[:]
        indx_to_extract += df_subset.cell.tolist()

    return sorted(indx_to_extract)


def get_coords(input_file, indx_to_extract):

    file_ext = os.path.basename(input_file).split('.')[-1]
    file_name = os.path.basename(input_file)[:-(len(file_ext)+1)]

    if file_ext == 'arc':
        with open(input_file, 'rt', encoding='latin-1') as file:
            lines = file.readlines()

        to_return = []
        for line in lines:
            if 'Empirical Formula' in line:
                natoms = int(line.split()[-2])
                break
        for i, line in enumerate(lines):
            if ('FINAL GEOMETRY OBTAINED' in line) and (int(lines[i+3].strip().split(':')[1]) in indx_to_extract):
                cell = int(lines[i+3].strip().split(':')[1])
                chunk = lines[i+4:i+4+natoms]
                sliced = []
                for c in chunk:
                    sliced.append(c.split())
                df = pd.DataFrame(sliced)
                # Creating a tuple (conf_name_Cell, coords)
                to_return.append((f"{file_name}_Cell_{cell}", df[[0, 1, 3, 5]]))
        return to_return

    elif file_ext == 'out':
        f = open(input_file, 'r')
        chunk = []
        cart = []

        # getting data from out
        # finding No. of atoms

        while True:
            line = f.readline()
            if "Empirical Formula" in line:
                natoms = int(line.split()[-2])
                f.close()
                break

        f = open(input_file, 'r')
        to_return = []
        while True:
            line = f.readline()
            if len(line) == 0:
                break

            if 79*'-' in line:
                while True:
                    line = f.readline()
                    if (79*'*' in line) or (len(line) == 0):
                        break

                    if 'CELL' in line:
                        cell = int(line.split(':')[1])

                    elif 'CARTESIAN COORDINATES' in line and cell in indx_to_extract:
                        utils.ignoreLines(f, 1)
                        cont = 0
                        chunk = []   
                        while cont < natoms:
                            chunk.append(f.readline())
                            cont += 1
                        cart = []
                        for c in chunk:
                            cart.append(c.split())
                            df = pd.DataFrame(cart)
                            to_return.append((f"{file_name}_Cell_{cell}", df[[1, 2, 3, 4]])) 
        return to_return
    else:
        raise ValueError(f"{input_file} does not have .arc or .out extension. Therefore is not readeable by ALEIMI.")


def main(
    input_file,
    boltzmann_file,
    energy_cut: float = 2,
    conformer_cut: float = 0,
    engine: str = 'psi4',
    machine: str = 'smaug',
    mkdir: bool = True,
    jobsh: bool = True,
    **keywords):

    if engine == 'psi4':
        InputExt = '.in'
    elif engine == 'orca':
        InputExt = '.inp'
    elif engine == 'gaussian':
        InputExt = '.gjf'
    else:
        print("Warning!: It was used 'in' as generic extension for the input "
              f"file for the non recognized engine: {engine}")
        InputExt = '.in'

    indx_to_extract = extract(boltzmann_file, energy_cut=energy_cut,
                              conformer_cut=conformer_cut)
    names_coords = get_coords(input_file, indx_to_extract)
    for name, coords in names_coords:
        INPUT_obj = templates.INPUT(engine, machine=machine, name=name, coords=coords, **keywords)
        if mkdir:
            utils.makedirs(name)
            INPUT_obj.write(os.path.join(name, f"{name}{InputExt}"), 'input')
            if jobsh:
                INPUT_obj.write(os.path.join(name, "job.sh"), 'jobsh')
        else:
            INPUT_obj.write(f"{name}{InputExt}", 'input')
            if jobsh:
                INPUT_obj.write(f"{name}.sh", 'jobsh')


if __name__ == '__main__':
    pass