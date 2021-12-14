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

if len(glob.glob('*.arc')) == 0:
    exit('No se puede encontrar el archivo .arc')


# =============================================================================
#      USER SPECIFICATIONS
# =============================================================================
#d_E = 0.001                                                         # eV
d_rmsd = 1.0                                                        # Angstrom
# =============================================================================

# reading files                                                               #
arcfs = glob.glob('*.arc')



def arc_reader (arcf):
    with open(arcf, 'rt', encoding='latin-1') as a:
        lines = a.readlines()
    # getting data from arc                                                   #
    totals_ev = []
    cells = []
    Class_E = []
    # finding No. of atoms
    for line in lines:
        if 'Empirical Formula' in line:
            natoms = int(line.split()[-2])
            break
    
    #CONTAINER = []
    CONTAINER__H = []
    #cart = []
    cart__H = []
    #atoms = []
    for i, line in enumerate(lines):

        if 'TOTAL ENERGY' in line:
            totals_ev.append(float(line.split()[-2]))
       
        elif 'Empirical Formula' in line:
            try:
                Class_E.append(float(lines[i+3].split('=')[1].split()[0]))
            except:
                Class_E.append('----')
            cells.append(int(lines[i+4].split(':')[1]))
    
        elif ('FINAL GEOMETRY OBTAINED' in line):
            chunk = lines[i+4:i+4+natoms]
            cart__H = []
            for c in chunk:
                #atoms.append(c.split()[0])
                #cart.append([c.split()[0].strip(), float(c.split()[1]), float(c.split()[3]), float(c.split()[5])])
                if c.split()[0] != 'H':
                    cart__H.append([c.split()[1], c.split()[3], c.split()[5]]) # No estoy tomando los atomos, solamente coordenadas c.split()[0].strip(),
            #CONTAINER.append(np.asarray(pd.DataFrame(cart)))
            CONTAINER__H.append(np.array(cart__H, dtype=np.float))
    #atoms = (np.asarray(atoms))
    # .... organizing
    paired = list(zip(cells, totals_ev, CONTAINER__H, Class_E)) # Esto genera un arreglo de tuplas, me une los arreglos
    ORDERED = sorted(paired, key=lambda x: x[1])  #Esto ordena la tupla segun la energia de menor a mayor
    return ORDERED #, atoms]


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
Ve minimos locales diferentes con energia parecida , pero no ve el mismo minimo reportado con un error o dos minimos diferentes en enrgia por el corte empleado pero diferentes en \
RMSD. Por ejemplo, una estructura que un cambio muy pequeño genera un cambio elevado en energia (dado por el cutoff empleado de enrgia que puede ser muy bajo).
Par DM quizas sea conveniente emplear esta version, pues no se necesita de tanta especificidad, y cambio muy pequeños no tendran una repesrcusion en los parametros bscados

En esta version solo evalua los RMSDs, quedamos con aquella molecula mas favorable energeticamente
pues de todas formas voy a a realizar calculos
de mayor nivel DFT
En el programa que me dieron hay algunos errores. el metodo de calculo de RMSD 
que usan no necesita una previa traslacion porque este lo tien implicito
"""
# =============================================================================

for arc in arcfs:
    name = arc.split('.')[0]
    ordered = (arc_reader(arc))
    


    for i, x in enumerate(range(len(ordered))):
        to_trash_degenerated = []
        for idx, y in enumerate(range(len(ordered))):
            if i < idx:
              
# =============================================================
#     CHECKING Geometric degeneracy
# =============================================================
                P = ordered[i][2]
                Q = ordered[idx][2]

                RMSD = rmsd.kabsch_rmsd(P, Q, translate = True)

                if RMSD <= d_rmsd:
                    # reject identical structure and kept the lowest energy (because ordered() is ordered using the enrgy, so idx always will have a grater enrgy)
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
