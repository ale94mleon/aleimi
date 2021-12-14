# -*- coding: utf-8 -*-
"""
===============================================================================
Created on    : Mon Aug 19 15:19:24 2019
Author        : Alejandro Martínez León
Mail          : [amleon@instec.cu, ale94mleon@gmail.com]
Affiliation   : Chemical Systems Modeling Group,
Affiliation   : Faculty of Radiochemistry, InSTEC-University of Havana, Cuba.
===============================================================================
DESCRIPTION   : Realiza el proceso de optimizacion clasico
DEPENDENCIES  :
===============================================================================
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import glob as glob
import pandas as pd
"""
Toma todos los archivos .mol, genera la cant de conformeros seleccionados
e imprime tanto los .sdf con todos los conformeros como .mop para optimizar en MOPAC
con el nivel de teoria seleccionado.
realiza opt con UFF
"""
if len(glob.glob('*.mol')) == 0:
    exit('No se puede encontrar el archivo .mol')
# =============================================================================
#      USER SPECIFICATIONS
# =============================================================================                                                         
d_rmsd = 0.2                                                        
Conf = 100
numThreads= 0
suppl = sorted(glob.glob('*.mol'))
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
def filter (mol, ids, opt, d_rmsd):
    print('Calculando RMSD ...')
    rmsd_reject=[]
    opt_reject =[]
    rmsd_no_reject=[]
    no_opt=0
    ids=list(ids)
    for i in ids:
        
        if opt[i][0]==1:
            no_opt +=1
            opt_reject.append(ids[i])
        else:
            for idx in ids:
                if i < idx:
                    RMSD = AllChem.GetConformerRMS(mol, i, idx, prealigned=True)
                    if RMSD <= d_rmsd:
                        #print(i,opt[i][1], idx, opt[idx][1])
                       #print(RMSD)
                        if opt[i][1] < opt[idx][1]:
                            rmsd_reject.append(ids[i])
                            rmsd_no_reject.append(ids[idx])
                        else:
                            rmsd_reject.append(ids[idx])
                            rmsd_no_reject.append(ids[i])
                            
    first_reject = opt_reject + rmsd_reject
    reject = no_repeat(first_reject)
    rmsd_no_reject = no_repeat(rmsd_no_reject)
    rmsd_reject = no_repeat(rmsd_reject)
    if len(rmsd_reject) !=0:
        print('Los confórmeros con Id %s no se imprimirán en el .mop por ser iguales (Según el error de RMSD = %s) a %s respectivamente y menos favorables energéticamente que los últimos mencionados.' % (rmsd_reject,d_rmsd,rmsd_no_reject))
    if len(opt_reject) !=0:
        print('Los confórmeros con Id %s no se imprimirán en el .mop pues no lograron ser optimizados..' % (opt_reject))

    to_use = []
    for item in ids:
        if item not in reject:
            to_use.append(item)
            
    return [to_use, no_opt, len(reject)-no_opt]

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

# =============================================================================
#    Comprobamos que los SMILES esten bien
# expuestas con la funcion filter y exportamos .sdf
# =============================================================================
Errors_mols = []
for item in suppl:
    try:
        posibble_molecule = (Chem.MolFromMolFile(item))
        posibble_molecule.GetNumAtoms()
    except:
        Errors_mols.append(item.split('.')[0])
if len(Errors_mols) != 0:
    exit('ERROR: Las moléculas : %s no son estructuras .mol correctas para RDKit. Por favor, retirelas de la carpeta de trabajo.' % (str(Errors_mols)))


molecules = []
for molecule in suppl:  
    molecules.append(Chem.MolFromMolFile(molecule))

# =============================================================================
#      Generamos conformaciones que cumplen con las condiciones de rmsd
# expuestas con la funcion filter y exportamos .sdf
# =============================================================================
print('Resumen\n')

picture=[]
for i, molecule in enumerate(molecules):
    
    print('Se están generarando %d conformaciones para la molécula %s ...' % (Conf, suppl[i].split('.')[0]))
    
    molecule = Chem.AddHs(molecule)
    natoms = molecule.GetNumAtoms()
    ids = AllChem.EmbedMultipleConfs(molecule, numConfs=Conf, numThreads=numThreads)
    opt = AllChem.UFFOptimizeMoleculeConfs(molecule, numThreads=numThreads, ignoreInterfragInteractions=False, maxIters=500)
    writer = Chem.SDWriter(suppl[i].split('.')[0]+'.sdf')
    to_use = filter(molecule, ids, opt, d_rmsd)
    
    for id in to_use[0]:
        writer.write(molecule, confId=id)
    writer.close()
    print(suppl[i].split('.')[0]+':\n  Error en opt: %d\n  Mismo RMSD: %d' % (to_use[1], to_use[2]))
    
# =============================================================================
#     Generando los .mop de los .sdf
# =============================================================================                                                         
    sdf = suppl[i].split('.')[0]+'.sdf'
                                                       
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
            comments = sdf.split('.')[0] + '   E = ' + str(round(opt[cell[cont]][1], 3)) + '  Molécula en el opt: '+str(cell[cont])                 
            final.write(comments+'\n')                 
            final.write('CELL: %d\n' % (cont+1))
            to_print.to_string(final, header=False, index=False)
            final.write('\n0\n')
    print('  Número de átomos: %d\n'% (natoms))
    final.close()
    
# =============================================================================
#exportar las imagenes de la mejor manera posible. en columnas de 1,2,3,4 o 5
# =============================================================================                                                         

ms = [x for x in molecules if x is not None]
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

img=Draw.MolsToGridImage(ms,molsPerRow=best,subImgSize=(700,700), legends=[ x.split('.')[0] for x in suppl])
img.save('Molecules used.png')
