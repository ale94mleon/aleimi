from aleimi import confgen, boltzmann, extractor, tools
suppl = "lactate.mol"

confgen.main(suppl, numConfs = 10, d_rmsd = 0.2, optimization = True, numThreads = 0, mopac_kewords =  'PM7 precise ef xyz geo-ok t=3h EPS=78.4')
tools.mopac('conf_lactate.mop')
boltzmann.main('conf_lactate.out', 1, out_path = 'conf_lactate.boltzmann')
extractor.main('conf_lactate.out', 'conf_lactate.boltzmann')

