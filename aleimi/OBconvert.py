# -*- coding: utf-8 -*-
"""
===============================================================================
Created on    : Sun Aug 25 18:49:11 2019
Author        : Alejandro Martínez León
Mail          : [amleon@instec.cu, ale94mleon@gmail.com]
Affiliation   : Chemical Systems Modeling Group,
Affiliation   : Faculty of Radiochemistry, InSTEC-University of Havana, Cuba.
===============================================================================
DESCRIPTION   : Usar con moderación. OpenBabel en dependencia del input pude dar discimiles errores.
DEPENDENCIES  : openbabel, glob
===============================================================================
"""


from openbabel import openbabel as ob
import glob as glob
import os, tempfile
from aleimi import tools



def obconvert(inpath, outpath, AddHydrogens = False, obabel = 'obabel'):
    """Convert  molecule ussing openbabel

    Args:
        input (str, path): input molecule.
        output (str, path): must have the extention of the molecule.
    """
    in_ext = os.path.basename(inpath).split('.')[-1]
    out_ext = os.path.basename(outpath).split('.')[-1] 

    if in_ext == 'pdbqt':
        tmpfile = tempfile.NamedTemporaryFile(suffix='.pdb')
        if AddHydrogens:
            tools.run(f"""
                    {obabel} {inpath} -O {tmpfile.name} -h > /dev/null 2>&1
                    {obabel} {tmpfile.name} -O {outpath} -h > /dev/null 2>&1
                    """)
        else:
            tools.run(f"""
                    {obabel} {inpath} -O {tmpfile.name} > /dev/null 2>&1
                    {obabel} {tmpfile.name} -O {outpath} > /dev/null 2>&1
                    """)
        tools.rm(tmpfile.name)


    else:
        obConversion = ob.OBConversion()
        obConversion.SetInAndOutFormats(in_ext, out_ext)
        mol = ob.OBMol()
        obConversion.ReadFile(mol, inpath)   # Open Babel will uncompressed automatically
        if AddHydrogens:
            mol.AddHydrogens()
        obConversion.WriteFile(mol, outpath)


if __name__ == '__main__':

    suppl = glob.glob('*')
    for i, item in enumerate(suppl):
        if '.py' in item:
            del suppl[i]
    for item in suppl:    
        obconvert(item, 'pdb')
