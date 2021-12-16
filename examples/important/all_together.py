from aleimi import confgen, boltzmann, extractor, tools

suppl = "lactate.mol"
confgen.main(suppl, numConfs = 100, rdkit_d_RMSD = 0.2, UFF = True, rdkit_numThreads = 0, mopac_keywords =  'PM7 precise ef xyz geo-ok t=3h EPS=78.4')
tools.mopac('conf_lactate.mop')
boltzmann.main('conf_lactate.out', 1, BOutPath = True)
extractor.main('conf_lactate.out', 'conf_lactate.boltzmann', energy_cut=1)

