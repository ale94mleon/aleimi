# -*- coding: utf-8 -*-
"""
===============================================================================
Created on    : Thu Aug 22 22:23:48 2019
Author        : Alejandro Martínez León
Mail          : [amleon@instec.cu, ale94mleon@gmail.com]
Affiliation   : Chemical Systems Modeling Group,
Affiliation   : Faculty of Radiochemistry, InSTEC-University of Havana, Cuba.
===============================================================================
DESCRIPTION   :
DEPENDENCIES  :
===============================================================================
"""
import glob
import rmsd
import numpy as np
import pandas as pd

if len(glob.glob('*.out')) == 0:
    exit('No se puede encontrar el archivo .out')


# =============================================================================
#      USER SPECIFICATIONS
# =============================================================================
d_E = 0.001                                                         # eV
d_rmsd = 1.0                                                        # Angstrom
# =============================================================================

# reading files                                                               #
outs = glob.glob('*.out')

def ignoreLines(f, n):
    for i in range(n): f.readline()

def out_reader (out):
    f = open(out, 'r')
    chunk = []
    totals_ev = []
    cells = []
    Class_E = []
    CONTAINER__H = []
    cart__H = []


    # getting data from out                                                       #
    # finding No. of atoms
    while True:
        l = f.readline()
        if "Empirical Formula" in l:
            natoms = int(l.split()[-2])
            break
            f.close

        
    f = open(out, 'r')
    while True:
        l = f.readline()
        k = l
        if len(l) == 0:break
               
        if '-------------------------------------------------------------------------------' in l:
            while True:
                k = f.readline()
                if ('*******************************************************************************' in k) or (len(k) == 0):break
                if 'CELL' in k: 
                    cells.append(int(k.split(':')[1]))
                    Class_E.append('----')
                elif 'TOTAL ENERGY' in k:
                    totals_ev.append(float(k.split()[-2]))
                elif 'CARTESIAN COORDINATES' in k:
                    ignoreLines(f, 1)
                    cont = 0
                    chunk = []        
                    while cont < natoms:
                        chunk.append(f.readline())
                        cont += 1
                    cart__H = []
                    for c in chunk:
                        if c.split()[1] != 'H':
                            cart__H.append([c.split()[2], c.split()[3], c.split()[4]]) # No estoy tomando los atomos, solamente coordenadas c.split()[0].strip(),

            CONTAINER__H.append(np.array(cart__H, dtype=np.float))
    f.close
    # .... organizing
    paired = list(zip(cells, totals_ev, CONTAINER__H, Class_E)) # Esto genera un arreglo de tuplas, me une los arreglos
    ORDERED = sorted(paired, key=lambda x: x[1])  #Esto ordena la tupla segun la energia de menor a mayor
    return ORDERED
    



# =============================================================================
#      CHECKING Energy degeneracy
"""
En el programa anterior la idea es:
    se parecene en enrgia:
        Si:
            Se parecen en geometria:
                Si: descarto una moleula ; No me quedo con ambas
        No: No hago nada, las moleculas son diferentes
Este lo que pretende es solamente evaluar los parametros geometricos puesto
que a mi entender el codigo anterior no se si tiene en cuenta el error del semiempirico,
pues dos estructuras que se parecen energeticamente pueden ser diferentes geometricamnete.
Me parece que solo ve minimos locales diferentes con energia parecida , pero no ve el mismo minimo reportado con un error.


En esta version solo evalua los RMSDs, quedamos con aquella molecula mas favorable energeticamente
pues de todas formas voy a a realizar calculos
de mayor nivel DFT
En el programa que me dieron hay algunos errores. el metodo de calculo de RMSD 
que usan no necesita una previa traslacion porque este lo tien implicito
"""
# =============================================================================

for out in outs:
    name = out.split('.')[0]
    ordered = (out_reader(out))
    


    for i, x in enumerate(range(len(ordered))):
        to_trash_degenerated = []
        for idx, y in enumerate(range(len(ordered))):
            if i < idx:
                Ei = ordered[i][1]
                Eidx = ordered[idx][1]
                delta = abs(Ei - Eidx)
                if delta <= d_E:
# =============================================================
#     CHECKING Geometric degeneracy
# =============================================================
                    P = ordered[i][2]
                    Q = ordered[idx][2]

                    RMSD = rmsd.kabsch_rmsd(P, Q, translate = True)

                    if RMSD <= d_rmsd:
                    # reject identical structure
                        to_trash_degenerated.append(idx)
# =========================================================================
#     FOR EACH STRUCTURE, eliminate degenerated and save lot of time
# =========================================================================
        to_trash_degenerated = sorted(to_trash_degenerated, reverse=True)
        [ordered.pop(x) for x in to_trash_degenerated]
        


# =============================================================================
#      WORKING with UNDEGENERATED. Cambie la manera de calculos los parametros:
#Me base en: James B. Foresman - Exploring Chemistry With Electronic Structure Methods 3rd edition (2015) pag 182
# y Mortimer_Physical Chemistry_(3rd.ed.-2008) pag 1045    
# =============================================================================
    k = 0.0019872                              # Boltzmann Constant (kcal/mol*K)
    T = 298.15                                 # Absolute T (K)
    DF = pd.DataFrame(columns=['cell','Class_E', 'ev', 'kcal', 'Emin_Ei',
                               'qi__Pi/Pmin__e^(Emin_Ei)/kT','Fraction_%__100*qi/q'])
    cells_und = [ordered[i][0] for i, x in enumerate(ordered)]
    Class_E = [ordered[i][3] for i, x in enumerate(ordered)]
    energ_ev = [ordered[i][1] for i, x in enumerate(ordered)]
    energ_kcal = [x * 23 for x in energ_ev]
    min_energ_kcal = min(energ_kcal)
    relative_kcal = [min_energ_kcal - x for x in energ_kcal]
    qi = [np.e**(E_r/(k*T)) for E_r in relative_kcal]
    q = sum(qi)
    Fraction = [100*i/q for i in qi]
    #Z = [np.e**(-(E/(k*T))) for E in energ_kcal] #no pudo calcular Z: verflowError: (34, 'Result too large') 
    #Pi_b = [(np.e**-(E/(k*T)))/Z for E in energ_kcal]
    # =============================================================================
    #     DATAFRAME
    # =============================================================================
    DF['cell'] = cells_und
    DF['Class_E'] = Class_E
    DF['ev'] = energ_ev
    DF['kcal'] = energ_kcal
    DF['Emin_Ei'] = relative_kcal
    DF['qi__Pi/Pmin__e^(Emin_Ei)/kT'] = qi
    DF['Fraction_%__100*qi/q'] = Fraction
    #DF['Pi_b = (e^((-Ei)/kT))/Z'] = Pi_b

    
    with open(name +'.boltzmann', 'wt') as rp:
        DF.to_string(rp)