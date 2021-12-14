# -*- coding: utf-8 -*-
"""
===============================================================================
Created on    : Mon Aug 19 15:19:24 2019
Author        : Alejandro Martínez León
Mail          : [amleon@instec.cu, ale94mleon@gmail.com]
Affiliation   : Chemical Systems Modeling Group,
Affiliation   : Faculty of Radiochemistry, InSTEC-University of Havana, Cuba.
===============================================================================
DESCRIPTION   :
DEPENDENCIES  :
===============================================================================
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import glob as glob
import pandas as pd
"""
Toma de un archivo .smi los SMILES, genera la cant de conformeros seleccionados
e imprime tanto los .sdf con todos los conformeros como .mop para optimizar en MOPAC
con el nivel de teoria seleccionado.
"""
if len(glob.glob('*.smi')) == 0:
    exit('No se puede encontrar el archivo .smi')
# =============================================================================
#      USER SPECIFICATIONS
# =============================================================================                                                         
d_rmsd = 0.2                                                        
Conf = 10
numThreads= 0
suppl = glob.glob('*.smi')[0]
keywords = 'PM7 precise ef xyz geo-ok t=3h EPS=78.4'
# =============================================================================

# =============================================================================
#      CHECKING geometry degeneracy
# La primera iteracion que corre por i busca que este optimizada la estructura
# La segunda busca que no se geometricamente la misma estructura
#Pueda suceder que alguna idx cumpla que no es geometricamente la misma estructura 
#y no la ve el segundo ciclo
#pero el primero cuando llegue a ella la vera y tomara su id
# =============================================================================
def filter (mol, ids, d_rmsd):
    print('Calculando RMSD ...')
    ides=list(ids)
    print(i)
    for i in ides:
        
        to_trash_degenerated = []
        for idx in ides:
            if i < idx:
                RMSD = AllChem.GetConformerRMS(mol, i, idx, prealigned=True)
                if RMSD <= d_rmsd:
                    to_trash_degenerated.append(idx)
        
        to_trash_degenerated = sorted(to_trash_degenerated, reverse=True)# se tiene que oredenar de mayor a mento para poder eliminar los elemnetos y que la indexacioin no cambie
        [ides.pop(x) for x in to_trash_degenerated]
           
            
    return [ides, len(ids)-len(ides)]

# =============================================================================
#      Esta funcion elimina los elementos repetidos de una lista y la ordena
# =============================================================================                                                         

def no_repeat(list):
    result = []
    for item in list:
        if item not in result:
            result.append(item)
            result = sorted(result)
    return result    

	
	

with open(suppl, 'rt') as file:
    smiles = file.readlines()
	
# =============================================================================
#    Comprobamos que los SMILES esten bien
# expuestas con la funcion filter y exportamos .sdf
# =============================================================================
Errors_mols = []
for i, item in enumerate(smiles):
    try:
        posibble_molecule = (Chem.MolFromSmiles(item))
        posibble_molecule.GetNumAtoms()
    except:
        Errors_mols.append(i+1)
if len(Errors_mols) !=  0:
    exit('ERROR: Las moléculas con Ids: %s no son estructuras SMILES correctas para RDKit. Por favor, retirelas del .smi.' % (str(Errors_mols)))


print('Resumen\n')
picture=[]

# =============================================================================
#      Generamos conformaciones que cumplen con las condiciones de rmsd
# expuestas con la funcion filter y exportamos .sdf
# =============================================================================
for i, smile in enumerate(smiles):
    
    picture.append(Chem.MolFromSmiles(smile))
    
    print('Se están generarando %d conformaciones para la molécula con Id = %s ...' % (Conf, i+1))
    
    molecule = Chem.AddHs(Chem.MolFromSmiles(smile))
    natoms = molecule.GetNumAtoms()
    ids = AllChem.EmbedMultipleConfs(molecule, numConfs=Conf, numThreads=numThreads)
    writer = Chem.SDWriter('conf_mol_'+str(i+1)+'.sdf')
    to_use = filter(molecule, ids, d_rmsd)
    
    for id in to_use[0]:
        writer.write(molecule, confId=id)
    writer.close()
    print('conf_mol_'+str(i+1)+':\n  Mismo RMSD: %d' % (to_use[1]))
    
# =============================================================================
#     Generando los .mop
# =============================================================================                                                         
    sdf = 'conf_mol_'+str(i+1)+'.sdf'
                                                       
    print('  Archivos de salida: %s, %s.mop' % (sdf, sdf.split('.')[0]))

    with open(sdf, 'rt') as file:
        lines = file.readlines()
    
    final=open(sdf.split('.')[0]+'.mop','w')
    
    cell = to_use[0] #Toma los valores reales que se usaron para imprimir el .sdf y asi buscar los opt de verdad, pq si se elimino una estructura tengo que tenerlo en cuenta
    cont = -1 #Lo inicializao en -1 pq la cantidad de veces que se encontrara la palabra RDKit coincide con la cantidad con el len(to_use[0])
    #todo esto para que en el .mop me saque el valor de energia optenido en la optimizacion
    for k, line in enumerate(lines):        
        if 'RDKit          3D' in line:
            cont += 1
            chunk = lines[(k+3):(k+3+(natoms))] 
            sliced = []
            for c in chunk:
                sliced.append(c.split())
            df = pd.DataFrame(sliced)
            to_print = df[[3, 0, 1, 2]]
            final.write(keywords+'\n')
            comments = sdf.split('.')[0] + '  Molécula en el sdf: '+str(cell[cont])                 
            final.write(comments+'\n')                 
            final.write('CELL: %d\n' % (cont+1))
            to_print.to_string(final, header=False, index=False)
            final.write('\n0\n')
    print('  Número de átomos: %d\n'% (natoms))
    final.close()
    
# =============================================================================
#exportar las imagenes de la mejor manera posible. en columnas de 1,2,3,4 o 5
# =============================================================================                                                         

ms = [x for x in picture if x is not None]
for m in ms: tmp=AllChem.Compute2DCoords(m)

cant = len(ms)
if cant==1:
    best=1
else:
    div = range(2,6)    
    tup = [[i - cant%i,i] for i in div]
    
    for t in tup:
        if t[0] == t[1]:
            t[0] = 0
    
    minimum_rest = min(tup, key=lambda x: x[0])
    for t in tup:
        if t[0] == minimum_rest[0] and t[1]>minimum_rest[1]:
            minimum_rest = t
    best = minimum_rest[1]



img=Draw.MolsToGridImage(ms,molsPerRow=best,subImgSize=(700,700), legends=[ 'Mol: '+str(x+1) for x in range (len(ms))])
img.save('Molecules used.png')
