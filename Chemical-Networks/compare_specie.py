import pandas as pd

first_network="xdr-new/species.dat"

second_network="complete_xdr/species.dat"

firstSpeciesList=list(pd.read_csv(first_network,header=None)[1].str.upper())
secondSpeciesList=list(pd.read_csv(second_network,header=None)[1].str.upper())

print([species for species in firstSpeciesList if species not in secondSpeciesList])
print([species for species in secondSpeciesList if species not in firstSpeciesList])