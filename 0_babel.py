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
# =============================================================================
#      USER SPECIFICATIONS
# =============================================================================                                                         
out_ext = 'pdb'
# =============================================================================




import openbabel as ob
import glob as glob


suppl = glob.glob('*')
for i, item in enumerate(suppl):
    if '.py' in item:
        del suppl[i]

def convert(name_file, out_ext):
    name = name_file.split('.')[0]
    in_ext = name_file.split('.')[1].lower()
    out_ext = out_ext.lower()
    print(name)
    print(in_ext)
    print(out_ext)

    obConversion = ob.OBConversion()
    obConversion.SetInAndOutFormats(in_ext, out_ext)
    mol = ob.OBMol()
    obConversion.ReadFile(mol, name_file)   # Open Babel will uncompress automatically
    #mol.AddHydrogens()
    obConversion.WriteFile(mol, name + '.' + out_ext)

for item in suppl:    
    convert(item, out_ext)
