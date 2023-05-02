from ase.phasediagram import PhaseDiagram

refs=[()] # add your computed strcutures and there related formation energies (eV/molecular unit) here
ref_E=[('Li',0),('C',0),('N',0)]
refs.extend(ref_E)
pd = PhaseDiagram(refs)
pd.plot(show=True,dims=2)
