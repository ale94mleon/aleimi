#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from aleimi import confgen, boltzmann, extractor, utils
import os, yaml, tempfile

cwd = os.getcwd()
tmp_path = tempfile.TemporaryDirectory()
os.chdir(tmp_path.name)
suppl1 = 'suppl1.smi'
suppl2 = 'suppl2.smi'

with open(suppl1, 'w') as s:
    s.write("COC(=O)C=1C=CC(=CC1)S(=O)(=O)N\nCCO")
with open(suppl2, 'w') as s:
    s.write("COC(=O)C\nCCO\nCNCO")

def test_module():

    mol_names = confgen.main('suppl1.smi')
    for mol_name in mol_names:
        print(mol_name)
        utils.mopac(f"{mol_name}.mop")
        boltzmann.main(f"{mol_name}.arc")
        extractor.main(f"{mol_name}.arc",f"{mol_name}_boltzmann.csv")

def test_cli_default_kwargs():
    utils.run('aleimi-run suppl2.smi')

def test_cli_user_kwargs():
    # Save the config as a yaml file
    params = {
    # aleimi.confgen.main
    'numConfs': 5,
    'rdkit_d_RMSD': 0.1,
    'UFF': True,
    'rdkit_numThreads': 0,
    'mopac_keywords': 'PM6 precise ef xyz geo-ok t=3h',

    # aleimi.boltzmann.main
    'Bd_rmsd': 2.0,
    'Bd_E': 0.0001,
    'BOutPath': True,
    
    # aleimi.extractor.main
    "energy_cut": 3,
    "conformer_cut": 4,
    "engine": 'psi4',
    }
    with open('params.yml', 'w') as p:
        yaml.dump(params, p)
    utils.run('aleimi-run suppl2.smi --params params.yml')

if __name__ == '__main__':
    test_module()