#! /usr/bin/python

import math
import os
import string
import struct
import sys
import time

# MakeRates version identifier
versionString = '2013-12-17'

# Use the Tkinter and Ttk modules, if available, to create a GUI interface
try:
    from tkinter import *
    import tkFileDialog
    import tkFont
    useTkinter = True
    try:
        from tkinter.ttk import *
        useTtk = True
    except:
        useTtk = False
except:
    try:
        from Tkinter import *
        import tkFileDialog
        import tkFont
        useTkinter = True
    except:
        useTkinter = False
    try:
        from ttk import *
        useTtk = True
    except:
        useTtk = False

# Default options
reactionFile = ''
speciesFile  = ''
outputPrefix = ''
sortSpecies = True
logForm = False
fileFormat = 'Rate05'
codeFormat = 'C'

# Reactant/product names to exclude from the species list
ignoreList = ['','#','e-']
ignoreList.extend(['ELECTR','PHOTON','CRP','CRPHOT']) # Gas-phase reactions
ignoreList.extend(['XRAY','XRSEC','XRLYA','XRPHOT'])  # X-ray reactions
ignoreList.extend(['FREEZE','CRH','PHOTD','THERM'])   # Gas-grain reactions
ignoreList.extend(['DESCR1','DESCR2','DEUVCR','DESOH2']) # Gas-grain reactions

# Determine the appropriate file format and read the entries in the specified reaction file
def read_reaction_file(fileName):
    input = open(fileName, mode='r')
    input.seek(0)
    data = input.readline().strip('\r\n')

    # Ignore header entries
    while len(data) < 5:
        data = input.readline().strip('\r\n')

    # Determine the format of the reaction file
    data = data.split(':')
    if len(data) > 10: fileFormat = 'rate12'
    else:
        data = data[0]
        data = data.split(',')
        if len(data) > 10: fileFormat = 'rate05'
        else:
            data = data[0]
            try:
                if len(data) < 102:
                    data += ' ' * (102 - len(data))
                data = data[:102]
                formatCode = '4sx8sx8sx8sx8sx8sx4sx4sx8sxxx5sxx8s1s5s5s1s4s'
                test = struct.unpack(formatCode, data)
                fileFormat = 'rate99'
            except:
                try:
                    formatCode = 'x4sx7sx7sx7sx7sx3sx3s8sx5sx8s1s4s4s'
                    test = struct.unpack(formatCode, data)
                    fileFormat = 'rate95'
                except:
                    if useTkinter:
                        app.statusMessage('ERROR! Unrecognised file format in input reaction file.', replace=True, error=True)
                    else:
                        sys.exit('\nERROR! Unrecognised file format in input reaction file\n')
                    return

    reactants = [] ; products = [] ; alpha = [] ; beta = [] ; gamma = [] ; labels = []
    while len(data) != 0 and data != ['']:
        if fileFormat == 'rate12':
            reactants.append([convert_species(value.strip()) for value in data[2:4]])
            products.append([convert_species(value.strip())  for value in data[4:8]])
            alpha.append(float(data[9]))
            beta.append(float(data[10]))
            gamma.append(float(data[11]))

            # Handle entries where the reference string contains a ':'
            if (len(data[16]) > 0 and data[16][-1] != '"'
            and len(data[17]) > 0 and data[17][0]  != '"'):
                data[16] += ':' + data.pop(17)

            labels.append([data[14],data[12],data[13],data[15],data[16]])

            # Account for missing label data (assume the maximum temperature range of 10 to 41000 K)
            if len(labels[-1]) < 4: labels[-1] = ['','10','41000','','']

            # Add an empty third reactant entry (since three reactants are not allowed in the Rate12 file format)
            if len(reactants[-1]) < 3: reactants[-1].append('')

            # Handle entries for the grain surface formation of H2 where the grain (labelled #) is included in the list of
            # products but not in the list of reactants (since three reactants are not allowed in the Rate12 file format)
            if reactants[-1][0] == 'H' and reactants[-1][1] == 'H' and reactants[-1][2] == '' and products[-1][0] == 'H2' and products[-1][1] == '#':
                reactants[-1][2] = '#'

            # Add entries for duplicate reactions valid over different temperature ranges
            for i in range(int(data[8])-1):
            	j = 17 + 9*i
                reactants.append(reactants[-1])
                products.append(products[-1])
                alpha.append(float(data[j+1]))
                beta.append(float(data[j+2]))
                gamma.append(float(data[j+3]))

                # Handle entries where the reference string contains a ':'
                if (len(data[j+8]) > 0 and data[j+8][-1] != '"'
                and len(data[j+9]) > 0 and data[j+9][0]  != '"'):
                    data[j+8] += ':' + data.pop(j+9)

                labels.append([data[j+6],data[j+4],data[j+5],data[j+7],data[j+8]])

            data = input.readline().strip('\r\n').split(':')

        elif fileFormat == 'rate05':
            reactants.append([convert_species(value.strip()) for value in data[1:4]])
            products.append([convert_species(value.strip())  for value in data[4:8]])
            alpha.append(float(data[8]))
            beta.append(float(data[9]))
            gamma.append(float(data[10]))
            labels.append([value for value in data[11:]])
            data = input.readline().strip('\r\n').split(',')

            # Account for missing label data (assume the maximum temperature range of 10 to 41000 K)
            if len(labels[-1]) < 4: labels[-1] = ['','10','41000','','']

        elif fileFormat == 'rate99':
            if len(data) < 102:
                data += ' ' * (102 - len(data))
            data = data[:102]
            data = struct.unpack(formatCode, data)
            reactants.append([convert_species(value.strip()) for value in data[1:4]])
            products.append([convert_species(value.strip())  for value in data[4:8]])
            alpha.append(float(data[8]))
            beta.append(float(data[9]))
            gamma.append(float(data[10]))
            labels.append([value.strip() for value in data[11:]])
            data = input.readline().strip('\r\n')

            # Account for missing label data (assume the maximum temperature range of 10 to 41000 K)
            if len(labels[-1]) < 4: labels[-1] = ['','10','41000','','']

        elif fileFormat == 'rate95':
            data = data[:77]
            data = struct.unpack(formatCode, data)
            reactants.append([convert_species(value.strip()) for value in data[1:3]])
            products.append([convert_species(value.strip())  for value in data[3:7]])
            alpha.append(float(data[7]))
            beta.append(float(data[8]))
            gamma.append(float(data[9]))

            # Account for missing label data (assume the maximum temperature range of 10 to 41000 K)
            labels.append([data[10].strip(),'10','41000'])
            labels[-1].append(['A','B','C','D','E',''][['1','2','3','4','5','9'].index(data[11][3])])
            labels[-1].append(data[12].strip())
            labels[-1].append(data[11].strip())

            # Add an empty third reactant entry (since three reactants are not allowed in the Rate12 file format)
            if len(reactants[-1]) < 3: reactants[-1].append('')

            data = input.readline().strip('\r\n')

    nReactions = len(reactants)
    input.close()
    return nReactions, reactants, products, alpha, beta, gamma, labels

# Determine the appropriate file format and read the entries in the specified species file
def read_species_file(fileName):
    input = open(fileName, mode='r')
    input.seek(0)
    data = input.readline().strip('\r\n')

    # Ignore header entries
    while len(data) < 5:
        data = input.readline().strip('\r\n')

    # Determine the format of the species file
    data = data.split(',')
    if len(data) > 3: fileFormat = 'rate05'
    else:
        data = data[0]
        try:
            formatCode = 'xx3sx8sxxx8sxx5s'
            test = struct.unpack(formatCode, data)
            fileFormat = 'rate99'
        except:
            try:
                formatCode = 'xxx3sxx9sxx3sxx9sxx3sxx9sxx3sxx9sxx3sxx9s'
                test = struct.unpack(formatCode, data)
                fileFormat = 'rate95'
            except:
                if useTkinter:
                    app.statusMessage('ERROR! Unrecognised file format in input species file.', replace=True, error=True)
                else:
                    sys.exit('\nERROR! Unrecognised file format in input species file\n')
                return

    species = [] ; abundance = [] ; mass = []
    while len(data) != 0 and data != ['']:
        if fileFormat == 'rate05':
            if data[1].lower() != 'e-' and data[1].upper() != 'ELECTR':
                species.append(convert_species(data[1]))
                abundance.append(float(data[2]))
                mass.append(float(data[3]))
            data = input.readline().split(',')
        elif fileFormat == 'rate99':
            data = struct.unpack(formatCode, data)
            if data[1].strip().lower() != 'e-' and data[1].strip().upper() != 'ELECTR':
                species.append(convert_species(data[1].strip()))
                abundance.append(float(data[2]))
                mass.append(float(data[3]))
            data = input.readline().strip('\r\n')
        elif fileFormat == 'rate95':
            data = data[8:]
            while len(data) > 0:
                if data[:9].strip().lower() != 'e-' and data[:9].strip().upper() != 'ELECTR':
                    species.append(convert_species(data[:9].strip()))
                    abundance.append(float(0))
                    mass.append(float(0))
                    if len(data) < 16: break
                    data = data[16:]
            data = input.readline().strip('\r\n')
    nSpecies = len(species)
    input.close()
    return nSpecies, species, abundance, mass

# Find all the species involved in a list of reactions
def find_all_species(reactantList, productList):
    speciesList = []
    for reactionReactants in reactantList:
        for reactant in reactionReactants:
            if ignoreList.count(reactant) == 0:
                if speciesList.count(reactant) == 0: speciesList.append(reactant)
    for reactionProducts in productList:
        for product in reactionProducts:
            if ignoreList.count(product) == 0:
                if speciesList.count(product) == 0: speciesList.append(product)
    nSpecies = len(speciesList)
    return nSpecies, speciesList

# Remove entries from a list of reactions when none of the listed species are involved
def find_all_reactions(speciesList, reactants, products, alpha, beta, gamma, labels):
    n = 0
    while n < len(reactants):
        for reactant in reactants[n]:
            if ignoreList.count(reactant) == 0:
                if speciesList.count(reactant) == 0:
                    reactants.pop(n) ; products.pop(n) ; alpha.pop(n) ; beta.pop(n) ; gamma.pop(n) ; labels.pop(n)
                    n -= 1
                    break
        n += 1
    n = 0
    while n < len(reactants):
        for product in products[n]:
            if ignoreList.count(product) == 0:
                if speciesList.count(product) == 0:
                    reactants.pop(n) ; products.pop(n) ; alpha.pop(n) ; beta.pop(n) ; gamma.pop(n) ; labels.pop(n)
                    n -= 1
                    break
        n += 1
    nReactions = len(reactants)
    return nReactions, reactants, products, alpha, beta, gamma, labels

# Return a list of "orphan" species that are either never formed or never destroyed
def check_orphan_species(speciesList, reactants, products):
    nSpecies = len(speciesList)

    # Count the total number of formation and destruction routes for each species
    nFormation = [] ; nDestruction = []
    if useTtk:
        app.progressValue.set(0)
        Style().configure('Horizontal.TProgressbar', background='#dcdad5')
    for n, species in enumerate(speciesList):
        if useTtk:
            app.progressValue.set(80.0*float(n+1)/float(nSpecies))
            app.status.update_idletasks()
        nFormation.append(sum([reactionProducts.count(species) for reactionProducts in products]))
        nDestruction.append(sum([reactionReactants.count(species) for reactionReactants in reactants]))

    # Check for species that are never formed or destroyed and produce an error message
    orphanList = [] ; missingList = []
    for n, species in enumerate(speciesList):
        if useTtk:
            app.progressValue.set(20.0*float(n+1)/float(nSpecies)+80.0)
            app.status.update_idletasks()
        if nFormation[n] == 0 and nDestruction[n] == 0:
             missingList.append(species)
        elif nFormation[n] == 0:
            orphanList.append(species)
            if useTkinter:
                app.statusMessage('\n\nERROR! Species "'+species+'" has destruction reaction(s) but no formation route.', error=True)
                print '\nERROR! Species "'+species+'" has destruction reaction(s) but no formation route.'
            else:
                sys.exit('\nERROR! Species "'+species+'" has destruction reaction(s) but no formation route\n')
        elif nDestruction[n] == 0:
            orphanList.append(species)
            if useTkinter:
                app.statusMessage('\n\nERROR! Species "'+species+'" has formation reaction(s) but no destruction route.', error=True)
                print '\nERROR! Species "'+species+'" has formation reaction(s) but no destruction route.'
            else:
                sys.exit('\nERROR! Species "'+species+'" has formation reaction(s) but no destruction route\n')
    if len(orphanList) > 0:
        return nFormation, nDestruction, missingList
    if len(missingList) > 0:
        if useTkinter:
            app.statusMessage('\n\nWARNING! The following species are missing from the reaction network:\n\n'+string.join(missingList,', '), error=True)
            print '\nWARNING! The following species are missing from the reaction network:\n'+string.join(missingList,', ')
        else:
            print '\nWARNING! The following species are missing from the reaction network:\n'+string.join(missingList,', ')
    return nFormation, nDestruction, missingList

# Convert a species name to the appropriate mixed-case form
def convert_species(species):
    upperCaseList = ['E-','ELECTR','H','D','HE','LI','C','N','O','F','NA','MG','SI','P','S','CL','CA','FE']
    mixedCaseList = ['e-','e-','H','D','He','Li','C','N','O','F','Na','Mg','Si','P','S','Cl','Ca','Fe']
    species = species.upper()
    for n in range(len(upperCaseList)):
        species = species.replace(upperCaseList[n], mixedCaseList[n])
    if len(species) > 1 and species[0] == 'M' and species[1] != 'g':
        species = '#' + species[1:]
    return species

# Sort a list of species first by their total number of destruction
# reactions and then by their total number of formation reactions
def sort_species(speciesList, abundanceList, nFormation, nDestruction):
    zippedList = zip(nDestruction, nFormation, speciesList, abundanceList)
    zippedList.sort()
    nDestruction, nFormation, speciesList, abundanceList = zip(*zippedList)
    zippedList = zip(nFormation, nDestruction, speciesList, abundanceList)
    zippedList.sort()
    nFormation, nDestruction, speciesList, abundanceList = zip(*zippedList)
    return speciesList, abundanceList

# Find the elemental constituents and molecular mass of each species in the supplied list
def find_constituents(speciesList):
    elementList = ['PAH','He','Li','Na','Mg','Si','Cl','Ca','Fe','H','D','C','N','O','F','P','S','#','+','-']
    elementMass = [420.0,4.0,7.0,23.0,24.0,28.0,35.0,40.0,56.0,1.0,2.0,12.0,14.0,16.0,19.0,31.0,32.0,0,0,0]
    nElements = len(elementList)
    speciesConstituents = []
    speciesMass = []

    for species in speciesList:
        constituents = []
        for element in elementList:
            constituents.append(0)
            for n in range(species.count(element)):
                index = species.index(element)+len(element)
                if species[index:index+2].isdigit():
                    constituents[-1] += int(species[index:index+2])
                    species = species[:index-len(element)]+species[index+2:]
                elif species[index:index+1].isdigit():
                    constituents[-1] += int(species[index:index+1])
                    species = species[:index-len(element)]+species[index+1:]
                else:
                    constituents[-1] += 1
                    species = species[:index-len(element)]+species[index:]

        # Calculate the total molecular mass as the sum of the elemental masses of each constituent
        speciesMass.append(sum([float(constituents[i])*float(elementMass[i]) for i in range(nElements)]))

        # Sort the elements in the constituent list by their atomic mass
        zippedList = zip(elementMass, elementList, constituents)
        zippedList.sort()
        sortedMasses, sortedElements, constituents = zip(*zippedList)
        speciesConstituents.append(constituents)

    # Sort the list of elements by their atomic mass
    zippedList = zip(elementMass, elementList)
    zippedList.sort()
    sortedMasses, sortedElements = zip(*zippedList)
    return speciesMass, speciesConstituents, sortedElements

# Check the conservation of elemental constituents and net charge in each reaction
def check_conservation(reactants, products, speciesList, elementList, speciesConstituents):
    failStatus = False
    nReactions = len(reactants)
    nElements = len(elementList)
    if useTtk:
        app.progressValue.set(0)
        Style().configure('Horizontal.TProgressbar', background='#dcdad5')

    # For each reaction, count the total number of each elemental constituent
    # contained in the reactants and in the products and check that they match
    for n in range(nReactions):
        reactantConstituents = [0]*nElements
        for reactant in reactants[n]:
            # Add the constituents of each reactant to the total constituent counts
            if reactant in speciesList:
                reactantConstituents = [reactantConstituents[i] + speciesConstituents[speciesList.index(reactant)][i] for i in range(nElements)]
        # Account for electrons, which are not included in the list of elements
        reactantConstituents[elementList.index('-')] += reactants[n].count('e-')
        # Determine the net charge of all reactants
        reactantConstituents.append(reactantConstituents[elementList.index('+')] - reactantConstituents[elementList.index('-')])

        productConstituents  = [0]*nElements
        for product in products[n]:
            # Add the constituents of each product to the total constituent counts
            if product in speciesList:
                productConstituents = [productConstituents[i] + speciesConstituents[speciesList.index(product)][i] for i in range(nElements)]
        # Account for electrons, which are not included in the list of elements
        productConstituents[elementList.index('-')] += products[n].count('e-')
        # Determine the net charge of all products
        productConstituents.append(productConstituents[elementList.index('+')] - productConstituents[elementList.index('-')])

        # Check that all constituent elements (and net charge) are conserved in the reaction
        for i in range(nElements):
            if i == elementList.index('#'): continue
            if i == elementList.index('+'): continue
            if i == elementList.index('-'): continue

            # Produce an error message if an element is not conserved (ignore freeze-out and gas-grain reactions)
            if productConstituents[i] != reactantConstituents[i] and reactants[n][1] != 'FREEZE' and reactantConstituents[elementList.index('#')] == 0:
                failStatus = True
                if useTkinter:
                    app.statusMessage('\n\nERROR! The total amount of "'+elementList[i]+'" is not conserved in the reaction:\n'
                                    + string.join(reactants[n][:3-reactants[n].count('')],' + ')+' --> '+string.join(products[n][:4-products[n].count('')],' + '), error=True)
                    print '\nERROR! The total amount of "'+elementList[i]+'" is not conserved in the reaction:\n' \
                           + string.join(reactants[n][:3-reactants[n].count('')],' + ')+' --> '+string.join(products[n][:4-products[n].count('')],' + ')
                else:
                    sys.exit('\nERROR! The total amount of "'+elementList[i]+'" is not conserved in the reaction '+str(n)+':\n' \
                           + string.join(reactants[n][:3-reactants[n].count('')],' + ')+' --> '+string.join(products[n][:4-products[n].count('')],' + ')+'\n')

        # Produce an error message if the net charge is not conserved (ignore freeze-out reactions)
        if productConstituents[-1] != reactantConstituents[-1] and reactants[n][1] != 'FREEZE':
            failStatus = True
            if useTkinter:
                app.statusMessage('\n\nERROR! The net charge is not conserved in the reaction:\n'
                                + string.join(reactants[n][:3-reactants[n].count('')],' + ')+' --> '+string.join(products[n][:4-products[n].count('')],' + '), error=True)
                print '\nERROR! The net charge is not conserved in the reaction:\n' \
                       + string.join(reactants[n][:3-reactants[n].count('')],' + ')+' --> '+string.join(products[n][:4-products[n].count('')],' + ')
            else:
                sys.exit('\nERROR! The net charge is not conserved in the reaction:\n' \
                       + string.join(reactants[n][:3-reactants[n].count('')],' + ')+' --> '+string.join(products[n][:4-products[n].count('')],' + ')+'\n')

        # Update the progress bar
        if useTtk:
            app.progressValue.set(100.0*float(n+1)/float(nReactions))
            app.status.update_idletasks()
    return failStatus

# Write the species file in the desired format
def write_species(fileName, speciesList, abundanceList, massList, fileFormat='rate05'):
    nSpecies = len(speciesList)
    output = open(fileName, mode='w')

    # Specify the appropriate format code for the desired file format
    if fileFormat == 'rate05': formatCode = '%i,%s,%.2E,%.1F\n'
    if fileFormat == 'rate99': formatCode = '  %3i %-8s   %8.2e  %5.1f\n'
    if fileFormat == 'rate95': formatCode = ' Y(%-3i)=%-8s'

    # Write the species names, abundances and molecular masses (if known)
    for n in range(nSpecies):
        if fileFormat == 'rate05': output.write(formatCode % (n+1,speciesList[n],abundanceList[n],massList[n]))
        if fileFormat == 'rate99': output.write(formatCode % (n+1,speciesList[n],abundanceList[n],massList[n]))
        if fileFormat == 'rate95': output.write(formatCode % (n+1,speciesList[n]))
        if fileFormat == 'rate95' and (n+1)%5 == 0: output.write('\n')

    # Include the electron as the final species
    if fileFormat == 'rate05': output.write(formatCode % (nSpecies+1,'e-',0,0))
    if fileFormat == 'rate99': output.write(formatCode % (0,'e-',0,0))
    if fileFormat == 'rate95': output.write('\n')
    output.close()

# Write the reaction file in the desired format
def write_reactions(fileName, reactants, products, alpha, beta, gamma, labels, fileFormat='rate05'):
    nReactions = len(reactants)
    output = open(fileName, mode='w')

    # Specify the appropriate format code for the desired file format
    if fileFormat == 'rate05': formatCode = '%i,%s,%s,%5.2E,%.2F,%.1F,%s\n'
    if fileFormat == 'rate99': formatCode = '%4i %-8s %-8s %-8s %-8s %-8s %-4s %-4s %8.2e   %5.2f  %8.1f%1s%5s%5s%1s%4s\n'
    if fileFormat == 'rate95': formatCode = ' %4i %-7s %-7s %-7s %-7s %-3s %-3s%8.2e %5.2f %8.1f%1s%4s%4s\n'

    # Write the reactants, products, Arrhenius equation parameters and accuracy labels (if available)
    if fileFormat == 'rate95': output.write('C\nC\nC\n')
    for n in range(nReactions):
        if fileFormat == 'rate05': output.write(formatCode % (n+1,string.join(reactants[n],','),string.join(products[n],','),alpha[n],beta[n],gamma[n],string.join(labels[n][0:5],',')))
        if fileFormat == 'rate99': output.write(formatCode % (n+1,reactants[n][0],reactants[n][1],reactants[n][2],products[n][0],products[n][1],products[n][2],products[n][3],alpha[n],beta[n],gamma[n],labels[n][0],labels[n][1],labels[n][2],labels[n][3],labels[n][4]))
        if fileFormat == 'rate95': output.write(formatCode % (n+1,reactants[n][0],reactants[n][1],products[n][0],products[n][1],products[n][2],products[n][3],alpha[n],beta[n],gamma[n],labels[n][0],(labels[n][5] if len(labels[n])==6 else '  A'+['1','2','3','4','5','9'][['A','B','C','D','E',''].index(labels[n][3])]),labels[n][4]))
    output.close()

# Create the conservation term for the desired species (an element, electron or dust grain)
def conserve_species(species, speciesConstituents, codeFormat='C'):
    elementList = ['#','+','-','H','D','He','Li','C','N','O','F','Na','Mg','Si','P','S','Cl','Ca','Fe','PAH']
    nSpecies = len(speciesConstituents)
    conservationEquation = ''
    # Handle the special case of electrons (i.e., charge conservation with both anions and cations)
    if species == 'e-':
        indexPos = elementList.index('+')
        indexNeg = elementList.index('-')
        for n in range(nSpecies):
            if speciesConstituents[n][indexPos] > 0:
                if len(conservationEquation) > 0: conservationEquation += '+'
                if codeFormat == 'C':   conservationEquation += multiple(speciesConstituents[n][indexPos])+'x['+str(n)+']'
                if codeFormat == 'F90': conservationEquation += multiple(speciesConstituents[n][indexPos])+'Y('+str(n+1)+')'
                if codeFormat == 'F77': conservationEquation += multiple(speciesConstituents[n][indexPos])+'Y('+str(n+1)+')'
            if speciesConstituents[n][indexNeg] > 0:
                conservationEquation += '-'
                if codeFormat == 'C':   conservationEquation += multiple(speciesConstituents[n][indexNeg])+'x['+str(n)+']'
                if codeFormat == 'F90': conservationEquation += multiple(speciesConstituents[n][indexNeg])+'Y('+str(n+1)+')'
                if codeFormat == 'F77': conservationEquation += multiple(speciesConstituents[n][indexNeg])+'Y('+str(n+1)+')'
    else:
        index = elementList.index(species)
        for n in range(nSpecies):
            if speciesConstituents[n][index] > 0:
                if len(conservationEquation) > 0: conservationEquation += '+'
                if codeFormat == 'C':   conservationEquation += multiple(speciesConstituents[n][index])+'x['+str(n)+']'
                if codeFormat == 'F90': conservationEquation += multiple(speciesConstituents[n][index])+'Y('+str(n+1)+')'
                if codeFormat == 'F77': conservationEquation += multiple(speciesConstituents[n][index])+'Y('+str(n+1)+')'
    if len(conservationEquation) > 0:
        if codeFormat == 'C':   conservationEquation = '  x_e = '+conservationEquation+';\n'
        if codeFormat == 'F90': conservationEquation = '      Y('+str(nSpecies+1)+') = '+conservationEquation+'\n'
        if codeFormat == 'F77': conservationEquation = '      X(1)  = '+conservationEquation+'\n'
    else:
        if codeFormat == 'C':   conservationEquation = '  x_e = 0;\n'
        if codeFormat == 'F90': conservationEquation = '      Y('+str(nSpecies+1)+') = 0\n'
        if codeFormat == 'F77': conservationEquation = '      X(1)  = 0\n'
    if codeFormat == 'F77': conservationEquation = truncate_line(conservationEquation)
    return conservationEquation

# Create the equations for the additional parameters needed in certain X-ray reaction rates
def xray_parameters(speciesList, codeFormat='C'):
    # Find the index numbers for H, H2, He and e- in the species list
    indexH  = speciesList.index('H')
    indexH2 = speciesList.index('H2')
    indexHe = speciesList.index('He')
    indexEl = len(speciesList)

    # Create the code string to calculate the parameters zeta_H, zeta_H2 and zeta_He,
    # i.e., 1/(W_i.x_i) in equation D.12 of Meijerink & Spaans (2005, A&A, 436, 397)
    if codeFormat == 'C':   xrayParameterEquations = '\n  /* The X-ray secondary ionization rates depend on the mean energies\n   * required to ionize H or H2 in a neutral gas mixture, 1/(W_i*x_i) */\n  zeta_H  = 1.0/(39.8*(1.0+12.2*pow(x_e,0.866))*(x['+str(indexH)+']+1.89*x['+str(indexH2)+']));\n  zeta_H2 = 1.0/(41.9*(1.0+6.72*pow((1.83*x_e/(1.0+0.83*x_e)),0.824))*(x['+str(indexH2)+']+0.53*x['+str(indexH)+']));\n  zeta_He = 1.0/(487.*(1.0+12.5*pow(x_e,0.994))*(x['+str(indexHe)+']));\n'
    if codeFormat == 'F90': xrayParameterEquations = '\n!     The X-ray secondary ionization rates depend on the mean energies\n!     required to ionize H or H2 in a neutral gas mixture, 1/(W_i*x_i)\n      ZETA_H  = 1.0D0/(39.8D0*(1.0D0+12.2D0*Y('+str(indexEl+1)+')**0.866D0)*(Y('+str(indexH+1)+')+1.89D0*Y('+str(indexH2+1)+')))\n      ZETA_H2 = 1.0D0/(41.9D0*(1.0D0+6.72D0*(1.83D0*Y('+str(indexEl+1)+')/(1.0D0+0.83D0*Y('+str(indexEl+1)+')))**0.824D0)*(Y('+str(indexH2+1)+')+0.53D0*Y('+str(indexH+1)+')))\n      ZETA_HE = 1.0D0/(487.0D0*(1.0D0+12.5D0*Y('+str(indexEl+1)+')**0.994D0)*(Y('+str(indexHe+1)+')))\n'
    if codeFormat == 'F77': xrayParameterEquations = '\nC     The X-ray secondary ionization rates depend on the mean energies\nC     required to ionize H or H2 in a neutral gas mixture, 1/(W_i*x_i)\n      ZETA_H  = 1.0D0/(39.8D0*(1.0D0+12.2D0*X(1)**0.866D0)*(Y('+str(indexH+1)+')+1.89D0*Y('+str(indexH2+1)+')))\n      ZETA_H2 = 1.0D0/(41.9D0*(1.0D0+6.72D0*(1.83D0*X(1)/(1.0D0+0.83D0*X(1)))**0.824D0)*(Y('+str(indexH2+1)+')+0.53D0*Y('+str(indexH+1)+')))\n      ZETA_HE = 1.0D0/(487.0D0*(1.0D0+12.5D0*X(1)**0.994D0)*(Y('+str(indexHe+1)+')))\n'
    if codeFormat == 'F77': xrayParameterEquations = truncate_line(xrayParameterEquations)
    return xrayParameterEquations

# Determine if the specified reactants and products represent a gas-grain reaction
def is_gasgrain_reaction(reactants, products):
    # Check for reactions labelled by gas-grain mechanisms
    if reactants[1] in ['FREEZE','CRH','PHOTD','THERM']: return True
    if reactants[1] in ['DESCR1','DESCR2','DEUVCR','DESOH2']: return True

    # Check for H2 formation on grain surfaces
    nReactants = len([species for species in reactants if species != ''])
    nProducts  = len([species for species in products  if species != ''])
    if nReactants == 2 and nProducts == 1:
        if reactants[0] == 'H' and reactants[1] == 'H' and products[0] == 'H2': return True
    if nReactants == 3 and nProducts == 2:
        if reactants[0] == 'H' and reactants[1] == 'H' and reactants[2] == '#' and products[0] == 'H2' and products[1] == '#': return True

    return False

# Create the appropriate multiplication string for a given number
def multiple(number):
    if number == 1: return ''
    else: return str(number)+'*'

# Truncate long lines for use in fixed-format Fortran code
def truncate_line(input, codeFormat='F77', continuationCode=None):
    lineLength = 72
    result = ''
    while len(input) > lineLength:
        # Handle newline characters in the input string
        if input.rfind('\n',0,lineLength) != -1:
            index = input.rfind('\n',0,lineLength) + 1
            result += input[:index]
            input = input[index:]
            continue

        # Truncate the line at the last possible occurrence of '+', '-', '*', or '/'
        index = max([input.rfind('+',0,lineLength),input.rfind('-',0,lineLength),input.rfind('*',0,lineLength),input.rfind('/',0,lineLength)])
        if codeFormat == 'F90':
            if continuationCode != None: result += input[:index]+' '+continuationCode.strip()+'\n'
            else: result += input[:index]+' &\n'
        else:
            result += input[:index]+'\n'
        if continuationCode != None:
            input = continuationCode+input[index:]
        else:
            input = '     &       '+input[index:]
    result += input
    return result

# Write the ODEs file in C language format
def write_odes_c(fileName, speciesList, constituentList, reactants, products, logForm=False):
    nSpecies = len(speciesList)
    nReactions = len(reactants)
    output = open(fileName, mode='w')

    # Determine if X-ray reactions are present in the chemical network
    if sum([reactantList.count('XRAY')+reactantList.count('XRSEC') for reactantList in reactants]) > 0:
        xrayReactions = True
    else:
        xrayReactions = False

    # Find the index numbers for H, H2 and He in the species list
    if 'H'  in speciesList: indexH  = speciesList.index('H')
    if 'H2' in speciesList: indexH2 = speciesList.index('H2')
    if 'He' in speciesList: indexHe = speciesList.index('He')

    # Write the comments and function header
    if logForm:
        if xrayReactions:
            fileHeader = '/*=======================================================================\n\n User-supplied f (ODEs) routine. Compute the function ydot = f(t,y)\n\n-----------------------------------------------------------------------*/\n\n/* Header files with descriptions of the contents used */\n\n#include <math.h>                    /* Standard math functions */\n#include <cvode/cvode.h>             /* CVODE functions and constants */\n#include <cvode/cvode_dense.h>       /* Prototype for CVDense solver */\n#include <nvector/nvector_serial.h>  /* Serial N_Vector types, functions, macros */\n#include <sundials/sundials_dense.h> /* Definition of type DlsMat (dense matrix) */\n#include <sundials/sundials_types.h> /* Definition of type realtype */\n\n/*-----------------------------------------------------------------------*/\n\n/* Type definition for user-supplied data passed to the solver functions */\n\ntypedef struct {\n  realtype *rate, n_H, T_g, x_e;\n} *User_Data;\n\n/*-----------------------------------------------------------------------*/\n\nint f(realtype t, N_Vector y, N_Vector ydot, void *user_data)\n{\n  realtype x['+str(nSpecies)+'], *ode, *rate;\n  realtype n_H, x_e, loss, form, zeta_H, zeta_H2, zeta_He;\n  User_Data data;\n\n  /* Obtain the pointer to the ydot vector data array */\n  ode = NV_DATA_S(ydot);\n\n  /* Retrieve the array of reaction rate coefficients and\n   * the total number density from the user-supplied data */\n  data = (User_Data) user_data;\n  rate = data->rate;\n  n_H = data->n_H;\n\n  /* Convert the abundances from logarithmic to normal form */\n  for (int i = 0; i < neq; i++) {\n    x[i] = exp(NV_Ith_S(y,i));\n  }\n\n  /* The electron abundance is a conserved quantity, given by the sum\n   * of the abundances of all ionized species in the chemical network */\n'
        else:
            fileHeader = '/*=======================================================================\n\n User-supplied f (ODEs) routine. Compute the function ydot = f(t,y)\n\n-----------------------------------------------------------------------*/\n\n/* Header files with descriptions of the contents used */\n\n#include <math.h>                    /* Standard math functions */\n#include <cvode/cvode.h>             /* CVODE functions and constants */\n#include <cvode/cvode_dense.h>       /* Prototype for CVDense solver */\n#include <nvector/nvector_serial.h>  /* Serial N_Vector types, functions, macros */\n#include <sundials/sundials_dense.h> /* Definition of type DlsMat (dense matrix) */\n#include <sundials/sundials_types.h> /* Definition of type realtype */\n\n/*-----------------------------------------------------------------------*/\n\n/* Type definition for user-supplied data passed to the solver functions */\n\ntypedef struct {\n  realtype *rate, n_H, T_g, x_e;\n} *User_Data;\n\n/*-----------------------------------------------------------------------*/\n\nint f(realtype t, N_Vector y, N_Vector ydot, void *user_data)\n{\n  realtype x['+str(nSpecies)+'], *ode, *rate;\n  realtype n_H, x_e, loss, form;\n  User_Data data;\n\n  /* Obtain the pointer to the ydot vector data array */\n  ode = NV_DATA_S(ydot);\n\n  /* Retrieve the array of reaction rate coefficients and\n   * the total number density from the user-supplied data */\n  data = (User_Data) user_data;\n  rate = data->rate;\n  n_H = data->n_H;\n\n  /* Convert the abundances from logarithmic to normal form */\n  for (int i = 0; i < neq; i++) {\n    x[i] = exp(NV_Ith_S(y,i));\n  }\n\n  /* The electron abundance is a conserved quantity, given by the sum\n   * of the abundances of all ionized species in the chemical network */\n'
    else:
        if xrayReactions:
            fileHeader = '/*=======================================================================\n\n User-supplied f (ODEs) routine. Compute the function ydot = f(t,y)\n\n-----------------------------------------------------------------------*/\n\n/* Header files with descriptions of the contents used */\n\n#include <math.h>                    /* Standard math functions */\n#include <cvode/cvode.h>             /* CVODE functions and constants */\n#include <cvode/cvode_dense.h>       /* Prototype for CVDense solver */\n#include <nvector/nvector_serial.h>  /* Serial N_Vector types, functions, macros */\n#include <sundials/sundials_dense.h> /* Definition of type DlsMat (dense matrix) */\n#include <sundials/sundials_types.h> /* Definition of type realtype */\n\n/*-----------------------------------------------------------------------*/\n\n/* Type definition for user-supplied data passed to the solver functions */\n\ntypedef struct {\n  realtype *rate, n_H, T_g, x_e;\n} *User_Data;\n\n/*-----------------------------------------------------------------------*/\n\nint f(realtype t, N_Vector y, N_Vector ydot, void *user_data)\n{\n  realtype *x, *ode, *rate;\n  realtype n_H, x_e, loss, form, zeta_H, zeta_H2, zeta_He;\n  User_Data data;\n\n  /* Obtain pointers to the y and ydot vector data arrays */\n  x = NV_DATA_S(y);\n  ode = NV_DATA_S(ydot);\n\n  /* Retrieve the array of reaction rate coefficients and\n   * the total number density from the user-supplied data */\n  data = (User_Data) user_data;\n  rate = data->rate;\n  n_H = data->n_H;\n\n  /* The electron abundance is a conserved quantity, given by the sum\n   * of the abundances of all ionized species in the chemical network */\n'
        else:
            fileHeader = '/*=======================================================================\n\n User-supplied f (ODEs) routine. Compute the function ydot = f(t,y)\n\n-----------------------------------------------------------------------*/\n\n/* Header files with descriptions of the contents used */\n\n#include <math.h>                    /* Standard math functions */\n#include <cvode/cvode.h>             /* CVODE functions and constants */\n#include <cvode/cvode_dense.h>       /* Prototype for CVDense solver */\n#include <nvector/nvector_serial.h>  /* Serial N_Vector types, functions, macros */\n#include <sundials/sundials_dense.h> /* Definition of type DlsMat (dense matrix) */\n#include <sundials/sundials_types.h> /* Definition of type realtype */\n\n/*-----------------------------------------------------------------------*/\n\n/* Type definition for user-supplied data passed to the solver functions */\n\ntypedef struct {\n  realtype *rate, n_H, T_g, x_e;\n} *User_Data;\n\n/*-----------------------------------------------------------------------*/\n\nint f(realtype t, N_Vector y, N_Vector ydot, void *user_data)\n{\n  realtype *x, *ode, *rate;\n  realtype n_H, x_e, loss, form;\n  User_Data data;\n\n  /* Obtain pointers to the y and ydot vector data arrays */\n  x = NV_DATA_S(y);\n  ode = NV_DATA_S(ydot);\n\n  /* Retrieve the array of reaction rate coefficients and\n   * the total number density from the user-supplied data */\n  data = (User_Data) user_data;\n  rate = data->rate;\n  n_H = data->n_H;\n\n  /* The electron abundance is a conserved quantity, given by the sum\n   * of the abundances of all ionized species in the chemical network */\n'
    output.write(fileHeader)

    # Prepare and write the electron conservation equation
    output.write(conserve_species('e-', constituentList, codeFormat='C'))

    # If X-ray reactions are present, write the additional terms needed to calculate their rates
    if xrayReactions:
        output.write(xray_parameters(speciesList, codeFormat='C'))

    # Prepare and write the loss and formation terms for each ODE
    output.write('\n  /* The ODEs created by MakeRates begin here... */\n')
    if useTtk:
        app.progressValue.set(0)
        Style().configure('Horizontal.TProgressbar', background='#dcdad5')
    for n in range(nSpecies):
        if useTtk:
            app.progressValue.set(100.0*float(n)/float(nSpecies))
            app.status.update_idletasks()
        species = speciesList[n]
        lossString = '' ; formString = ''

        # Loss terms
        for i in range(nReactions):
            if reactants[i].count(species) > 0:
                lossString += '-'+multiple(reactants[i].count(species))+'rate['+str(i)+']'

                # Gas-grain reaction rates include a density factor to account for the number density of grains
                if is_gasgrain_reaction(reactants[i], products[i]):
                    lossString += '*n_H'
                    continue

                # Purely gas phase or grain surface reactions
                for reactant in speciesList:
                    if reactant == species:
                        for j in range(reactants[i].count(reactant)-1):
                            lossString += '*x['+str(speciesList.index(reactant))+']*n_H'
                    else:
                        for j in range(reactants[i].count(reactant)):
                            lossString += '*x['+str(speciesList.index(reactant))+']*n_H'

                # Multiply by the electron density when electrons appear as a reactant
                for j in range(reactants[i].count('e-')):
                    lossString += '*x_e*n_H'

                # X-ray induced secondary ionization
                if reactants[i].count('XRSEC') == 1:
                    if reactants[i].count('H') == 1:
                        lossString += '*zeta_H'
                    elif reactants[i].count('H2') == 1:
                        lossString += '*zeta_H2'
                    elif reactants[i].count('He') == 1:
                        lossString += '*zeta_He'
                    else:
                        lossString += '*zeta_H'

                # Photoreactions due to X-ray induced secondary photons (Lyman-alpha from excited H)
                if reactants[i].count('XRLYA') == 1:
                    lossString += '*x['+str(indexH)+']*zeta_H'

                # Photoreactions due to X-ray induced secondary photons (Lyman-Werner from excited H2)
                if reactants[i].count('XRPHOT') == 1:
                    lossString += '*x['+str(indexH2)+']*zeta_H2'

            # Formation terms
            if products[i].count(species) > 0:
                formString += '+'+multiple(products[i].count(species))+'rate['+str(i)+']'

                # Gas-grain reaction rates include a density factor to account for the number density of grains
                if is_gasgrain_reaction(reactants[i], products[i]):
                    formString += '*x['+str(speciesList.index(reactants[i][0]))+']*n_H'
                    continue

                # Purely gas phase or grain surface reactions
                for reactant in speciesList:
                    for j in range(reactants[i].count(reactant)):
                        formString += '*x['+str(speciesList.index(reactant))+']'

                # Multiply by the electron density when electrons appear as a reactant
                for j in range(reactants[i].count('e-')):
                    formString += '*x_e'

                if sum([speciesList.count(reactant) for reactant in reactants[i]]) > 1 or reactants[i].count('e-') > 0:
                    formString += '*n_H'

                # X-ray induced secondary ionization
                if reactants[i].count('XRSEC') == 1:
                    if reactants[i].count('H') == 1:
                        formString += '*zeta_H'
                    elif reactants[i].count('H2') == 1:
                        formString += '*zeta_H2'
                    elif reactants[i].count('He') == 1:
                        formString += '*zeta_He'
                    else:
                        formString += '*zeta_H'

                # Photoreactions due to X-ray induced secondary photons (Lyman-alpha from excited H)
                if reactants[i].count('XRLYA') == 1:
                    formString += '*x['+str(indexH)+']*zeta_H'

                # Photoreactions due to X-ray induced secondary photons (Lyman-Werner from excited H2)
                if reactants[i].count('XRPHOT') == 1:
                    formString += '*x['+str(indexH2)+']*zeta_H2'

        if lossString != '':
            lossString = '  loss = '+lossString+';\n'
            output.write(lossString)
        if formString != '':
            formString = '  form = '+formString+';\n'
            output.write(formString)
        ydotString = '  ode['+str(n)+'] = '
        if formString != '':
            ydotString += 'form'
            if lossString != '': ydotString += '+'
        if lossString != '':
            ydotString += 'x['+str(n)+']*loss'
        ydotString += ';\n'
        output.write(ydotString)

    # If the logarithmic form of the ODEs is to be used, divide each by its abundance
    if logForm:
        output.write('\n')
        output.write('\n  /* Convert the ODEs from dy/dt to d[ln(y)]/dt by dividing each by its abundance */\n')
        for n in range(nSpecies):
            output.write('  ode['+str(n)+'] = ode['+str(n)+']/x['+str(n)+'];\n')

    # Write the function footer
    fileFooter = '\n  /* Store the electron abundance in the user data */\n  data->x_e = x_e;\n\n  return(0);\n}\n/*=======================================================================*/\n'
    output.write(fileFooter)
    output.close()

# Write the ODEs file in F90 language format
def write_odes_f90(fileName, speciesList, constituentList, reactants, products):
    nSpecies = len(speciesList)
    nReactions = len(reactants)
    output = open(fileName, mode='w')

    # Determine if X-ray reactions are present in the chemical network
    if sum([reactantList.count('XRAY')+reactantList.count('XRSEC') for reactantList in reactants]) > 0:
        xrayReactions = True
    else:
        xrayReactions = False

    # Find the index numbers for H, H2 and He in the species list
    if 'H'  in speciesList: indexH  = speciesList.index('H')
    if 'H2' in speciesList: indexH2 = speciesList.index('H2')
    if 'He' in speciesList: indexHe = speciesList.index('He')

    # Write the comments and function header
    if xrayReactions:
        fileHeader = '      SUBROUTINE F(NEQ,T,Y,YDOT)\n\n      USE DEFINITIONS\n      USE HEALPIX_TYPES\n      USE MAINCODE_MODULE, ONLY : RATE,N_H\n\n      IMPLICIT NONE\n\n      INTEGER(KIND=I4B), INTENT(IN) :: NEQ\n      REAL(KIND=DP), INTENT(IN)     :: T\n      REAL(KIND=DP), INTENT(INOUT)  :: Y(1:NEQ+1)\n      REAL(KIND=DP), INTENT(OUT)    :: YDOT(1:NEQ)\n\n      REAL(KIND=DP) :: LOSS,PROD,ZETA_H,ZETA_H2,ZETA_HE\n\n'
    else:
        fileHeader = '      SUBROUTINE F(NEQ,T,Y,YDOT)\n\n      USE DEFINITIONS\n      USE HEALPIX_TYPES\n      USE MAINCODE_MODULE, ONLY : RATE,N_H\n\n      IMPLICIT NONE\n\n      INTEGER(KIND=I4B), INTENT(IN) :: NEQ\n      REAL(KIND=DP), INTENT(IN)     :: T\n      REAL(KIND=DP), INTENT(INOUT)  :: Y(1:NEQ+1)\n      REAL(KIND=DP), INTENT(OUT)    :: YDOT(1:NEQ)\n\n      REAL(KIND=DP) :: LOSS,PROD\n\n'
    output.write(fileHeader)

    # Prepare and write the electron conservation equation
    output.write(conserve_species('e-', constituentList, codeFormat='F90'))

    # If X-ray reactions are present, write the additional terms needed to calculate their rates
    if xrayReactions:
        output.write(xray_parameters(speciesList, codeFormat='F90'))

    # Prepare and write the loss and formation terms for each ODE
    output.write('\n!     The ODEs created by MakeRates begin here...\n')
    if useTtk:
        app.progressValue.set(0)
        Style().configure('Horizontal.TProgressbar', background='#dcdad5')
    for n in range(nSpecies):
        if useTtk:
            app.progressValue.set(100.0*float(n)/float(nSpecies))
            app.status.update_idletasks()
        species = speciesList[n]
        lossString = '' ; formString = ''

        # Loss terms
        for i in range(nReactions):
            if reactants[i].count(species) > 0:
                lossString += '-'+multiple(reactants[i].count(species))+'RATE('+str(i+1)+')'

                # Gas-grain reaction rates include a density factor to account for the number density of grains
                if is_gasgrain_reaction(reactants[i], products[i]):
                    lossString += '*N_H'
                    continue

                # Purely gas phase or grain surface reactions
                for reactant in speciesList:
                    if reactant == species:
                        for j in range(reactants[i].count(reactant)-1):
                            lossString += '*Y('+str(speciesList.index(reactant)+1)+')*N_H'
                    else:
                        for j in range(reactants[i].count(reactant)):
                            lossString += '*Y('+str(speciesList.index(reactant)+1)+')*N_H'

                # Multiply by the electron density when electrons appear as a reactant
                for j in range(reactants[i].count('e-')):
                    lossString += '*Y('+str(nSpecies+1)+')*N_H'

                # X-ray induced secondary ionization
                if reactants[i].count('XRSEC') == 1:
                    if reactants[i].count('H') == 1:
                        lossString += '*ZETA_H'
                    elif reactants[i].count('H2') == 1:
                        lossString += '*ZETA_H2'
                    elif reactants[i].count('He') == 1:
                        lossString += '*ZETA_HE'
                    else:
                        lossString += '*ZETA_H'

                # Photoreactions due to X-ray induced secondary photons (Lyman-alpha from excited H)
                if reactants[i].count('XRLYA') == 1:
                    lossString += '*Y('+str(indexH+1)+')*ZETA_H'

                # Photoreactions due to X-ray induced secondary photons (Lyman-Werner from excited H2)
                if reactants[i].count('XRPHOT') == 1:
                    lossString += '*Y('+str(indexH2+1)+')*ZETA_H2'

            # Formation terms
            if products[i].count(species) > 0:
                formString += '+'+multiple(products[i].count(species))+'RATE('+str(i+1)+')'

                # Gas-grain reaction rates include a density factor to account for the number density of grains
                if is_gasgrain_reaction(reactants[i], products[i]):
                    formString += '*Y('+str(speciesList.index(reactants[i][0])+1)+')*N_H'
                    continue

                # Purely gas phase or grain surface reactions
                for reactant in speciesList:
                    for j in range(reactants[i].count(reactant)):
                        formString += '*Y('+str(speciesList.index(reactant)+1)+')'

                # Multiply by the electron density when electrons appear as a reactant
                for j in range(reactants[i].count('e-')):
                    formString += '*Y('+str(nSpecies+1)+')'

                if sum([speciesList.count(reactant) for reactant in reactants[i]]) > 1 or reactants[i].count('e-') > 0:
                    formString += '*N_H'

                # X-ray induced secondary ionization
                if reactants[i].count('XRSEC') == 1:
                    if reactants[i].count('H') == 1:
                        formString += '*ZETA_H'
                    elif reactants[i].count('H2') == 1:
                        formString += '*ZETA_H2'
                    elif reactants[i].count('He') == 1:
                        formString += '*ZETA_HE'
                    else:
                        formString += '*ZETA_H'

                # Photoreactions due to X-ray induced secondary photons (Lyman-alpha from excited H)
                if reactants[i].count('XRLYA') == 1:
                    formString += '*Y('+str(indexH+1)+')*ZETA_H'

                # Photoreactions due to X-ray induced secondary photons (Lyman-Werner from excited H2)
                if reactants[i].count('XRPHOT') == 1:
                    formString += '*Y('+str(indexH2+1)+')*ZETA_H2'

        if lossString != '':
            lossString = '      LOSS = '+lossString+'\n'
            output.write(lossString)
        if formString != '':
            formString = '      PROD = '+formString+'\n'
            output.write(formString)
        ydotString = '      YDOT('+str(n+1)+') = '
        if formString != '':
            ydotString += 'PROD'
            if lossString != '': ydotString += '+'
        if lossString != '':
            ydotString += 'Y('+str(n+1)+')*LOSS'
        ydotString += '\n'
        output.write(ydotString)

    # Write the function footer
    fileFooter = '\n      RETURN\n      END\n'
    output.write(fileFooter)
    output.close()

# Write the ODEs file in F77 language format
def write_odes_f77(fileName, speciesList, constituentList, reactants, products):
    nSpecies = len(speciesList)
    nReactions = len(reactants)
    output = open(fileName, mode='w')

    # Determine if X-ray reactions are present in the chemical network
    if sum([reactantList.count('XRAY')+reactantList.count('XRSEC') for reactantList in reactants]) > 0:
        xrayReactions = True
    else:
        xrayReactions = False

    # Find the index numbers for H, H2 and He in the species list
    if 'H'  in speciesList: indexH  = speciesList.index('H')
    if 'H2' in speciesList: indexH2 = speciesList.index('H2')
    if 'He' in speciesList: indexHe = speciesList.index('He')

    # Write the comments and function header
    if xrayReactions:
        fileHeader = "      SUBROUTINE F(NEQ,T,Y,YDOT)\n      INCLUDE 'header.f'\n      INTEGER          NEQ\n      DOUBLE PRECISION T,Y,YDOT\n      DOUBLE PRECISION D,LOSS,PROD\n      DOUBLE PRECISION ZETA_H,ZETA_H2,ZETA_HE\n      DIMENSION        NEQ(*),Y(*),YDOT(*)\n\nC     Set D to the gas density for use in the ODEs\n      D=DENS(DEPTH)\n\n"
    else:
        fileHeader = "      SUBROUTINE F(NEQ,T,Y,YDOT)\n      INCLUDE 'header.f'\n      INTEGER          NEQ\n      DOUBLE PRECISION T,Y,YDOT\n      DOUBLE PRECISION D,LOSS,PROD\n      DIMENSION        NEQ(*),Y(*),YDOT(*)\n\nC     Set D to the gas density for use in the ODEs\n      D=DENS(DEPTH)\n\n"
    output.write(fileHeader)

    # Prepare and write the electron conservation equation
    output.write(conserve_species('e-', constituentList, codeFormat='F77'))

    # If X-ray reactions are present, write the additional terms needed to calculate their rates
    if xrayReactions:
        output.write(xray_parameters(speciesList, codeFormat='F77'))

    # Prepare and write the loss and formation terms for each ODE
    output.write('\nC     The ODEs created by MakeRates begin here...\n')
    if useTtk:
        app.progressValue.set(0)
        Style().configure('Horizontal.TProgressbar', background='#dcdad5')
    for n in range(nSpecies):
        if useTtk:
            app.progressValue.set(100.0*float(n)/float(nSpecies))
            app.status.update_idletasks()
        species = speciesList[n]
        lossString = '' ; formString = ''

        # Loss terms
        for i in range(nReactions):
            if reactants[i].count(species) > 0:
                lossString += '-'+multiple(reactants[i].count(species))+'K('+str(i+1)+')'

                # Gas-grain reaction rates include a density factor to account for the number density of grains
                if is_gasgrain_reaction(reactants[i], products[i]):
                    lossString += '*D'
                    continue

                # Purely gas phase or grain surface reactions
                for reactant in speciesList:
                    if reactant == species:
                        for j in range(reactants[i].count(reactant)-1):
                            lossString += '*Y('+str(speciesList.index(reactant)+1)+')*D'
                    else:
                        for j in range(reactants[i].count(reactant)):
                            lossString += '*Y('+str(speciesList.index(reactant)+1)+')*D'

                # Multiply by the electron density when electrons appear as a reactant
                for j in range(reactants[i].count('e-')):
                    lossString += '*X(1)*D'

                # X-ray induced secondary ionization
                if reactants[i].count('XRSEC') == 1:
                    if reactants[i].count('H') == 1:
                        lossString += '*ZETA_H'
                    elif reactants[i].count('H2') == 1:
                        lossString += '*ZETA_H2'
                    elif reactants[i].count('He') == 1:
                        lossString += '*ZETA_HE'
                    else:
                        lossString += '*ZETA_H'

                # Photoreactions due to X-ray induced secondary photons (Lyman-alpha from excited H)
                if reactants[i].count('XRLYA') == 1:
                    lossString += '*Y('+str(indexH+1)+')*ZETA_H'

                # Photoreactions due to X-ray induced secondary photons (Lyman-Werner from excited H2)
                if reactants[i].count('XRPHOT') == 1:
                    lossString += '*Y('+str(indexH2+1)+')*ZETA_H2'

            # Formation terms
            if products[i].count(species) > 0:
                formString += '+'+multiple(products[i].count(species))+'K('+str(i+1)+')'

                # Gas-grain reaction rates include a density factor to account for the number density of grains
                if is_gasgrain_reaction(reactants[i], products[i]):
                    formString += '*Y('+str(speciesList.index(reactants[i][0])+1)+')*D'
                    continue

                # Purely gas phase or grain surface reactions
                for reactant in speciesList:
                    for j in range(reactants[i].count(reactant)):
                        formString += '*Y('+str(speciesList.index(reactant)+1)+')'

                # Multiply by the electron density when electrons appear as a reactant
                for j in range(reactants[i].count('e-')):
                    formString += '*X(1)'

                if sum([speciesList.count(reactant) for reactant in reactants[i]]) > 1 or reactants[i].count('e-') > 0:
                    formString += '*D'

                # X-ray induced secondary ionization
                if reactants[i].count('XRSEC') == 1:
                    if reactants[i].count('H') == 1:
                        formString += '*ZETA_H'
                    elif reactants[i].count('H2') == 1:
                        formString += '*ZETA_H2'
                    elif reactants[i].count('He') == 1:
                        formString += '*ZETA_HE'
                    else:
                        formString += '*ZETA_H'

                # Photoreactions due to X-ray induced secondary photons (Lyman-alpha from excited H)
                if reactants[i].count('XRLYA') == 1:
                    formString += '*Y('+str(indexH+1)+')*ZETA_H'

                # Photoreactions due to X-ray induced secondary photons (Lyman-Werner from excited H2)
                if reactants[i].count('XRPHOT') == 1:
                    formString += '*Y('+str(indexH2+1)+')*ZETA_H2'

        if lossString != '':
            lossString = '      LOSS = '+lossString+'\n'
            lossString = truncate_line(lossString)
            output.write(lossString)
        if formString != '':
            formString = '      PROD = '+formString+'\n'
            formString = truncate_line(formString)
            output.write(formString)
        ydotString = '      YDOT('+str(n+1)+') = '
        if formString != '':
            ydotString += 'PROD'
            if lossString != '': ydotString += '+'
        if lossString != '':
            ydotString += 'Y('+str(n+1)+')*LOSS'
        ydotString += '\n'
        ydotString = truncate_line(ydotString)
        output.write(ydotString)

    # Write the function footer
    fileFooter = '\n      RETURN\n      END\n'
    output.write(fileFooter)
    output.close()

# Write the Jacobian matrix file in C language format
def write_jac_c(fileName, speciesList, reactants, products, logForm=False):
    nSpecies = len(speciesList)
    nReactions = len(reactants)
    output = open(fileName, mode='w')

    # Determine if X-ray reactions are present in the chemical network
    if sum([reactantList.count('XRAY')+reactantList.count('XRSEC') for reactantList in reactants]) > 0:
        xrayReactions = True
    else:
        xrayReactions = False

    # Find the index numbers for H, H2 and He in the species list
    if 'H'  in speciesList: indexH  = speciesList.index('H')
    if 'H2' in speciesList: indexH2 = speciesList.index('H2')
    if 'He' in speciesList: indexHe = speciesList.index('He')

    # Specify the appropriate format code to represent the matrix indices
    if nSpecies >= 1000:
        formatCode = '%4i'
    elif nSpecies >=100:
        formatCode = '%3i'
    elif nSpecies >= 10:
        formatCode = '%2i'
    else:
        formatCode = '%i'

    # Write the comments and function header
    if xrayReactions:
        fileHeader = '/*=======================================================================\n\n User-supplied Jacobian routine. Compute the function J(t,y) = df/dy\n\n-----------------------------------------------------------------------*/\n\n/* Header files with descriptions of the contents used */\n\n#include <math.h>                    /* Standard math functions */\n#include <cvode/cvode.h>             /* CVODE functions and constants */\n#include <cvode/cvode_dense.h>       /* Prototype for CVDense solver */\n#include <nvector/nvector_serial.h>  /* Serial N_Vector types, functions, macros */\n#include <sundials/sundials_dense.h> /* Definition of type DlsMat (dense matrix) */\n#include <sundials/sundials_types.h> /* Definition of type realtype */\n\n/*-----------------------------------------------------------------------*/\n\n/* Type definition for user-supplied data passed to the solver functions */\n\ntypedef struct {\n  realtype *rate, n_H, T_g, x_e;\n} *User_Data;\n\n/*-----------------------------------------------------------------------*/\n\nint Jac(long int N, realtype t, N_Vector y, N_Vector fy, DlsMat J,\n        void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)\n{\n  realtype *x, *rate;\n  realtype n_H, x_e, zeta_H, zeta_H2, zeta_He;\n  User_Data data;\n\n  /* Obtain a pointer to the y vector data array */\n  x = NV_DATA_S(y);\n\n  /* Retrieve the array of reaction rate coefficients, total number\n   * density and the electron abundance from the user-supplied data */\n  data = (User_Data) user_data;\n  rate = data->rate;\n  n_H = data->n_H;\n  x_e = data->x_e;\n'
    else:
        fileHeader = '/*=======================================================================\n\n User-supplied Jacobian routine. Compute the function J(t,y) = df/dy\n\n-----------------------------------------------------------------------*/\n\n/* Header files with descriptions of the contents used */\n\n#include <math.h>                    /* Standard math functions */\n#include <cvode/cvode.h>             /* CVODE functions and constants */\n#include <cvode/cvode_dense.h>       /* Prototype for CVDense solver */\n#include <nvector/nvector_serial.h>  /* Serial N_Vector types, functions, macros */\n#include <sundials/sundials_dense.h> /* Definition of type DlsMat (dense matrix) */\n#include <sundials/sundials_types.h> /* Definition of type realtype */\n\n/*-----------------------------------------------------------------------*/\n\n/* Type definition for user-supplied data passed to the solver functions */\n\ntypedef struct {\n  realtype *rate, n_H, T_g, x_e;\n} *User_Data;\n\n/*-----------------------------------------------------------------------*/\n\nint Jac(long int N, realtype t, N_Vector y, N_Vector fy, DlsMat J,\n        void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)\n{\n  realtype *x, *rate;\n  realtype n_H, x_e;\n  User_Data data;\n\n  /* Obtain a pointer to the y vector data array */\n  x = NV_DATA_S(y);\n\n  /* Retrieve the array of reaction rate coefficients, total number\n   * density and the electron abundance from the user-supplied data */\n  data = (User_Data) user_data;\n  rate = data->rate;\n  n_H = data->n_H;\n  x_e = data->x_e;\n'
    output.write(fileHeader)

    # If X-ray reactions are present, write the additional terms needed to calculate their rate partial derivatives
    if xrayReactions:
        output.write(xray_parameters(speciesList, codeFormat='C'))
        additionalString = '\n  /* The additional partial derivative terms for the X-ray secondary ionization reactions begin here... */\n'

    # Prepare and write the terms for each Jacobian matrix element
    output.write('\n  /* The Jacobian matrix created by MakeRates begin here... */\n')
    if useTtk:
        app.progressValue.set(0)
        Style().configure('Horizontal.TProgressbar', background='#dcdad5')
    for n in range(nSpecies):
        if useTtk:
            app.progressValue.set(100.0*float(n+1)/float(nSpecies))
            app.status.update_idletasks()
        species1 = speciesList[n]
        for m in range(nSpecies):
            species2 = speciesList[m]
            matrixString = ''

            # Loss terms for species1
            for i in range(nReactions):
                if reactants[i].count(species1) > 0 and reactants[i].count(species2) > 0:

                    # Gas-grain reaction rates include a density factor to account for the number density of grains
                    if is_gasgrain_reaction(reactants[i], products[i]):
                        matrixString += '-'+multiple(reactants[i].count(species1))+'rate['+str(i)+']*n_H'
                        continue

                    matrixString += '-'+multiple(reactants[i].count(species1))+multiple(reactants[i].count(species1))+'rate['+str(i)+']'

                    # Purely gas phase or grain surface reactions
                    for reactant in speciesList:
                        if reactant == species2:
                            for j in range(reactants[i].count(reactant)-1):
                                matrixString += '*x['+str(speciesList.index(reactant))+']*n_H'
                        else:
                            for j in range(reactants[i].count(reactant)):
                                matrixString += '*x['+str(speciesList.index(reactant))+']*n_H'

                    # Multiply by the electron density when electrons appear as a reactant
                    for j in range(reactants[i].count('e-')):
                        matrixString += '*x_e*n_H'

                    # X-ray induced secondary ionization
                    if reactants[i].count('XRSEC') == 1:
                        if reactants[i].count('H') == 1:
                            matrixString = matrixString[:-len('-rate['+str(i)+']')]
                            additionalString += '  J->cols['+str(formatCode % indexH2)+']['+str(formatCode % n)+'] -= rate['+str(i)+']*x['+str(indexH) +']*zeta_H*(-1.89/(x['+str(indexH)+']+1.89*x['+str(indexH2)+']));\n'
                            additionalString += '  J->cols['+str(formatCode % indexH) +']['+str(formatCode % n)+'] -= rate['+str(i)+']*x['+str(indexH2)+']*zeta_H*(+1.89/(x['+str(indexH)+']+1.89*x['+str(indexH2)+']));\n'
                        elif reactants[i].count('H2') == 1:
                            matrixString = matrixString[:-len('-rate['+str(i)+']')]
                            additionalString += '  J->cols['+str(formatCode % indexH2)+']['+str(formatCode % n)+'] -= rate['+str(i)+']*x['+str(indexH) +']*zeta_H2*(+0.53/(x['+str(indexH2)+']+0.53*x['+str(indexH)+']));\n'
                            additionalString += '  J->cols['+str(formatCode % indexH) +']['+str(formatCode % n)+'] -= rate['+str(i)+']*x['+str(indexH2)+']*zeta_H2*(-0.53/(x['+str(indexH2)+']+0.53*x['+str(indexH)+']));\n'
                        elif reactants[i].count('He') == 1:
                            matrixString += '*zeta_He'
                        else:
                            matrixString += '*zeta_H'
                            additionalString += '  J->cols['+str(formatCode % indexH2)+']['+str(formatCode % n)+'] -= rate['+str(i)+']*x['+str(n)+']*zeta_H*(-1.89/(x['+str(indexH)+']+1.89*x['+str(indexH2)+']));\n'
                            additionalString += '  J->cols['+str(formatCode % indexH) +']['+str(formatCode % n)+'] -= rate['+str(i)+']*x['+str(n)+']*zeta_H*(-1.00/(x['+str(indexH)+']+1.89*x['+str(indexH2)+']));\n'

                    # Photoreactions due to X-ray induced secondary photons (Lyman-alpha from excited H)
                    if reactants[i].count('XRLYA') == 1:
                        matrixString += '*x['+str(indexH)+']*zeta_H'
#                        additionalString += '  J->cols['+str(formatCode % indexH2)+']['+str(formatCode % n)+'] -= rate['+str(i)+']*x['+str(n)+']*zeta_H*(-1.89*x['+str(indexH)+ ']/(x['+str(indexH)+']+1.89*x['+str(indexH2)+']));\n'
#                        additionalString += '  J->cols['+str(formatCode % indexH) +']['+str(formatCode % n)+'] -= rate['+str(i)+']*x['+str(n)+']*zeta_H*(+1.89*x['+str(indexH2)+']/(x['+str(indexH)+']+1.89*x['+str(indexH2)+']));\n'

                    # Photoreactions due to X-ray induced secondary photons (Lyman-Werner from excited H2)
                    if reactants[i].count('XRPHOT') == 1:
                        matrixString += '*x['+str(indexH2)+']*zeta_H2'
#                        additionalString += '  J->cols['+str(formatCode % indexH2)+']['+str(formatCode % n)+'] -= rate['+str(i)+']*x['+str(n)+']*zeta_H2*(+0.53*x['+str(indexH)+ ']/(x['+str(indexH2)+']+0.53*x['+str(indexH)+']));\n'
#                        additionalString += '  J->cols['+str(formatCode % indexH) +']['+str(formatCode % n)+'] -= rate['+str(i)+']*x['+str(n)+']*zeta_H2*(-0.53*x['+str(indexH2)+']/(x['+str(indexH2)+']+0.53*x['+str(indexH)+']));\n'

            # Formation terms for species1
            for i in range(nReactions):
                if products[i].count(species1) > 0 and reactants[i].count(species2) > 0:
                    matrixString += '+'+multiple(products[i].count(species1))+'rate['+str(i)+']'

                    # Gas-grain reaction rates include a density factor to account for the number density of grains
                    if is_gasgrain_reaction(reactants[i], products[i]):
                        matrixString += '*n_H'
                        continue

                    # Purely gas phase or grain surface reactions
                    for reactant in speciesList:
                        if reactant == species2:
                            for j in range(reactants[i].count(reactant)-1):
                                matrixString += '*x['+str(speciesList.index(reactant))+']*n_H'
                        else:
                            for j in range(reactants[i].count(reactant)):
                                matrixString += '*x['+str(speciesList.index(reactant))+']*n_H'

                    # Multiply by the electron density when electrons appear as a reactant
                    for j in range(reactants[i].count('e-')):
                        matrixString += '*x_e*n_H'

                    # X-ray induced secondary ionization
                    if reactants[i].count('XRSEC') == 1:
                        if reactants[i].count('H') == 1:
                            matrixString = matrixString[:-len('+rate['+str(i)+']')]
                            additionalString += '  J->cols['+str(formatCode % indexH2)+']['+str(formatCode % n)+'] += rate['+str(i)+']*x['+str(indexH) +']*zeta_H*(-1.89/(x['+str(indexH)+']+1.89*x['+str(indexH2)+']));\n'
                            additionalString += '  J->cols['+str(formatCode % indexH) +']['+str(formatCode % n)+'] += rate['+str(i)+']*x['+str(indexH2)+']*zeta_H*(+1.89/(x['+str(indexH)+']+1.89*x['+str(indexH2)+']));\n'
                        elif reactants[i].count('H2') == 1:
                            matrixString = matrixString[:-len('+rate['+str(i)+']')]
                            additionalString += '  J->cols['+str(formatCode % indexH2)+']['+str(formatCode % n)+'] += rate['+str(i)+']*x['+str(indexH) +']*zeta_H2*(+0.53/(x['+str(indexH2)+']+0.53*x['+str(indexH)+']));\n'
                            additionalString += '  J->cols['+str(formatCode % indexH) +']['+str(formatCode % n)+'] += rate['+str(i)+']*x['+str(indexH2)+']*zeta_H2*(-0.53/(x['+str(indexH2)+']+0.53*x['+str(indexH)+']));\n'
                        elif reactants[i].count('He') == 1:
                            matrixString += '*zeta_He'
                        else:
                            matrixString += '*zeta_H'
                            additionalString += '  J->cols['+str(formatCode % indexH2)+']['+str(formatCode % n)+'] += rate['+str(i)+']*x['+str(m)+']*zeta_H*(-1.89/(x['+str(indexH)+']+1.89*x['+str(indexH2)+']));\n'
                            additionalString += '  J->cols['+str(formatCode % indexH) +']['+str(formatCode % n)+'] += rate['+str(i)+']*x['+str(m)+']*zeta_H*(-1.00/(x['+str(indexH)+']+1.89*x['+str(indexH2)+']));\n'

                    # Photoreactions due to X-ray induced secondary photons (Lyman-alpha from excited H)
                    if reactants[i].count('XRLYA') == 1:
                        matrixString += '*x['+str(indexH)+']*zeta_H'
#                        additionalString += '  J->cols['+str(formatCode % indexH2)+']['+str(formatCode % n)+'] += rate['+str(i)+']*x['+str(m)+']*zeta_H*(-1.89*x['+str(indexH)+ ']/(x['+str(indexH)+']+1.89*x['+str(indexH2)+']));\n'
#                        additionalString += '  J->cols['+str(formatCode % indexH) +']['+str(formatCode % n)+'] += rate['+str(i)+']*x['+str(m)+']*zeta_H*(+1.89*x['+str(indexH2)+']/(x['+str(indexH)+']+1.89*x['+str(indexH2)+']));\n'

                    # Photoreactions due to X-ray induced secondary photons (Lyman-Werner from excited H2)
                    if reactants[i].count('XRPHOT') == 1:
                        matrixString += '*x['+str(indexH2)+']*zeta_H2'
#                        additionalString += '  J->cols['+str(formatCode % indexH2)+']['+str(formatCode % n)+'] += rate['+str(i)+']*x['+str(m)+']*zeta_H2*(+0.53*x['+str(indexH)+ ']/(x['+str(indexH2)+']+0.53*x['+str(indexH)+']));\n'
#                        additionalString += '  J->cols['+str(formatCode % indexH) +']['+str(formatCode % n)+'] += rate['+str(i)+']*x['+str(m)+']*zeta_H2*(-0.53*x['+str(indexH2)+']/(x['+str(indexH2)+']+0.53*x['+str(indexH)+']));\n'

            if matrixString != '':
                matrixString = '  J->cols['+str(formatCode % m)+']['+str(formatCode % n)+'] = '+matrixString+';\n'
                output.write(matrixString)

    # If X-ray reactions are present, write their additional partial derivative terms
    if xrayReactions:
        output.write(additionalString)

    # If the logarithmic form of the ODEs is to be used, multiply each matrix element J_ij by y_j/y_i
    if logForm:
        output.write('\n')
        for n in range(nSpecies):
            for m in range(nSpecies):
                if n != m:
                    output.write('  J->cols['+str(formatCode % m)+']['+str(formatCode % n)+'] = J->cols['+str(m)+']['+str(n)+']*x['+str(m)+']/x['+str(n)+'];\n')

    # Write the function footer
    fileFooter = '\n  return(0);\n}\n/*=======================================================================*/\n'
    output.write(fileFooter)
    output.close()

# Write the Jacobian matrix file in F90 language format
def write_jac_f90(fileName, speciesList, reactants, products):
    nSpecies = len(speciesList)
    nReactions = len(reactants)
    output = open(fileName, mode='w')

    # Determine if X-ray reactions are present in the chemical network
    if sum([reactantList.count('XRAY')+reactantList.count('XRSEC') for reactantList in reactants]) > 0:
        xrayReactions = True
    else:
        xrayReactions = False

    # Find the index numbers for H, H2 and He in the species list
    if 'H'  in speciesList: indexH  = speciesList.index('H')
    if 'H2' in speciesList: indexH2 = speciesList.index('H2')
    if 'He' in speciesList: indexHe = speciesList.index('He')

    # Specify the appropriate format code to represent the matrix indices
    if nSpecies >= 1000:
        formatCode = '%4i'
    elif nSpecies >=100:
        formatCode = '%3i'
    elif nSpecies >= 10:
        formatCode = '%2i'
    else:
        formatCode = '%i'

    # Write the comments and function header
    if xrayReactions:
        fileHeader = '      SUBROUTINE JAC(NEQ,T,Y,ML,MU,PD,NROWPD)\n\n      USE DEFINITIONS\n      USE HEALPIX_TYPES\n      USE MAINCODE_MODULE, ONLY : N_H,RATE\n\n      IMPLICIT NONE\n\n      INTEGER(KIND=I4B), INTENT(IN) :: NEQ,ML,MU,NROWPD\n      REAL(KIND=DP), INTENT(IN)     :: T\n      REAL(KIND=DP), INTENT(IN)     :: Y(1:NEQ+1)\n      REAL(KIND=DP), INTENT(OUT)    :: PD(1:NROWPD,1:NEQ)\n\n      REAL(KIND=DP) :: ZETA_H,ZETA_H2,ZETA_HE\n'
    else:
        fileHeader = '      SUBROUTINE JAC(NEQ,T,Y,ML,MU,PD,NROWPD)\n\n      USE DEFINITIONS\n      USE HEALPIX_TYPES\n      USE MAINCODE_MODULE, ONLY : N_H,RATE\n\n      IMPLICIT NONE\n\n      INTEGER(KIND=I4B), INTENT(IN) :: NEQ,ML,MU,NROWPD\n      REAL(KIND=DP), INTENT(IN)     :: T\n      REAL(KIND=DP), INTENT(IN)     :: Y(1:NEQ+1)\n      REAL(KIND=DP), INTENT(OUT)    :: PD(1:NROWPD,1:NEQ)\n'
    output.write(fileHeader)

    # If X-ray reactions are present, write the additional terms needed to calculate their rate partial derivatives
    if xrayReactions:
        output.write(xray_parameters(speciesList, codeFormat='F90'))
        additionalString = '\n!     The additional partial derivative terms for the X-ray secondary ionization reactions begin here...\n'

    # Prepare and write the terms for each Jacobian matrix element
    output.write('\n!     The Jacobian matrix created by MakeRates begins here...\n')
    if useTtk:
        app.progressValue.set(0)
        Style().configure('Horizontal.TProgressbar', background='#dcdad5')
    for n in range(nSpecies):
        if useTtk:
            app.progressValue.set(100.0*float(n)/float(nSpecies))
            app.status.update_idletasks()
        species1 = speciesList[n]
        for m in range(nSpecies):
            species2 = speciesList[m]
            matrixString = ''

            # Loss terms for species1
            for i in range(nReactions):
                if reactants[i].count(species1) > 0 and reactants[i].count(species2) > 0:

                    # Gas-grain reaction rates include a density factor to account for the number density of grains
                    if is_gasgrain_reaction(reactants[i], products[i]):
                        matrixString += '-'+multiple(reactants[i].count(species1))+'RATE('+str(i+1)+')*N_H'
                        continue

                    matrixString += '-'+multiple(reactants[i].count(species1))+multiple(reactants[i].count(species1))+'RATE('+str(i+1)+')'

                    # Purely gas phase or grain surface reactions
                    for reactant in speciesList:
                        if reactant == species2:
                            for j in range(reactants[i].count(reactant)-1):
                                matrixString += '*Y('+str(speciesList.index(reactant)+1)+')*N_H'
                        else:
                            for j in range(reactants[i].count(reactant)):
                                matrixString += '*Y('+str(speciesList.index(reactant)+1)+')*N_H'

                    # Multiply by the electron density when electrons appear as a reactant
                    for j in range(reactants[i].count('e-')):
                        matrixString += '*Y('+str(nSpecies+1)+')*N_H'

                    # X-ray induced secondary ionization
                    if reactants[i].count('XRSEC') == 1:
                        if reactants[i].count('H') == 1:
                            matrixString = matrixString[:-len('-RATE('+str(i+1)+')')]
                            additionalString += '      PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH2+1))+') = PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH2+1))+')-RATE('+str(i+1)+')*Y('+str(indexH+1) +')*ZETA_H*(-1.89/(Y('+str(indexH+1)+')+1.89*Y('+str(indexH2+1)+')))\n'
                            additionalString += '      PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH+1)) +') = PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH+1)) +')-RATE('+str(i+1)+')*Y('+str(indexH2+1)+')*ZETA_H*(+1.89/(Y('+str(indexH+1)+')+1.89*Y('+str(indexH2+1)+')))\n'
                        elif reactants[i].count('H2') == 1:
                            matrixString = matrixString[:-len('-RATE('+str(i+1)+')')]
                            additionalString += '      PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH2+1))+') = PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH2+1))+')-RATE('+str(i+1)+')*Y('+str(indexH+1) +')*ZETA_H2*(+0.53/(Y('+str(indexH2+1)+')+0.53*Y('+str(indexH+1)+')))\n'
                            additionalString += '      PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH+1)) +') = PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH+1)) +')-RATE('+str(i+1)+')*Y('+str(indexH2+1)+')*ZETA_H2*(-0.53/(Y('+str(indexH2+1)+')+0.53*Y('+str(indexH+1)+')))\n'
                        elif reactants[i].count('He') == 1:
                            matrixString += '*ZETA_HE'
                        else:
                            matrixString += '*ZETA_H'
                            additionalString += '      PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH2+1))+') = PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH2+1))+')-RATE('+str(i+1)+')*Y('+str(n+1)+')*ZETA_H*(-1.89/(Y('+str(indexH+1)+')+1.89*Y('+str(indexH2+1)+')))\n'
                            additionalString += '      PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH+1)) +') = PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH+1)) +')-RATE('+str(i+1)+')*Y('+str(n+1)+')*ZETA_H*(-1.00/(Y('+str(indexH+1)+')+1.89*Y('+str(indexH2+1)+')))\n'

                    # Photoreactions due to X-ray induced secondary photons (Lyman-alpha from excited H)
                    if reactants[i].count('XRLYA') == 1:
                        matrixString += '*Y('+str(indexH+1)+')*ZETA_H'
#                        additionalString += TO BE ADDED!!!
#                        additionalString += TO BE ADDED!!!

                    # Photoreactions due to X-ray induced secondary photons (Lyman-Werner from excited H2)
                    if reactants[i].count('XRPHOT') == 1:
                        matrixString += '*Y('+str(indexH2+1)+')*ZETA_H2'
#                        additionalString += TO BE ADDED!!!
#                        additionalString += TO BE ADDED!!!

            # Formation terms for species1
            for i in range(nReactions):
                if products[i].count(species1) > 0 and reactants[i].count(species2) > 0:
                    matrixString += '+'+multiple(products[i].count(species1))+'RATE('+str(i+1)+')'

                    # Gas-grain reaction rates include a density factor to account for the number density of grains
                    if is_gasgrain_reaction(reactants[i], products[i]):
                        matrixString += '*N_H'
                        continue

                    # Purely gas phase or grain surface reactions
                    for reactant in speciesList:
                        if reactant == species2:
                            for j in range(reactants[i].count(reactant)-1):
                                matrixString += '*Y('+str(speciesList.index(reactant)+1)+')*N_H'
                        else:
                            for j in range(reactants[i].count(reactant)):
                                matrixString += '*Y('+str(speciesList.index(reactant)+1)+')*N_H'

                    # Multiply by the electron density when electrons appear as a reactant
                    for j in range(reactants[i].count('e-')):
                        matrixString += '*Y('+str(nSpecies+1)+')*N_H'

                    # X-ray induced secondary ionization
                    if reactants[i].count('XRSEC') == 1:
                        if reactants[i].count('H') == 1:
                            matrixString = matrixString[:-len('+RATE('+str(i+1)+')')]
                            additionalString += '      PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH2+1))+') = PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH2+1))+')+RATE('+str(i+1)+')*Y('+str(indexH+1) +')*ZETA_H*(-1.89/(Y('+str(indexH+1)+')+1.89*Y('+str(indexH2+1)+')))\n'
                            additionalString += '      PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH+1)) +') = PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH+1)) +')+RATE('+str(i+1)+')*Y('+str(indexH2+1)+')*ZETA_H*(+1.89/(Y('+str(indexH+1)+')+1.89*Y('+str(indexH2+1)+')))\n'
                        elif reactants[i].count('H2') == 1:
                            matrixString = matrixString[:-len('+rate['+str(i)+']')]
                            additionalString += '      PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH2+1))+') = PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH2+1))+')+RATE('+str(i+1)+')*Y('+str(indexH+1) +')*ZETA_H2*(+0.53/(Y('+str(indexH2+1)+')+0.53*Y('+str(indexH+1)+')))\n'
                            additionalString += '      PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH+1)) +') = PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH+1)) +')+RATE('+str(i+1)+')*Y('+str(indexH2+1)+')*ZETA_H2*(-0.53/(Y('+str(indexH2+1)+')+0.53*Y('+str(indexH+1)+')))\n'
                        elif reactants[i].count('He') == 1:
                            matrixString += '*ZETA_HE'
                        else:
                            matrixString += '*ZETA_H'
                            additionalString += '      PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH2+1))+') = PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH2+1))+')+RATE('+str(i+1)+')*Y('+str(m+1)+']*ZETA_H*(-1.89/(Y('+str(indexH+1)+')+1.89*Y('+str(indexH2+1)+')))\n'
                            additionalString += '      PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH+1)) +') = PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH+1)) +')+RATE('+str(i+1)+')*Y('+str(m+1)+']*ZETA_H*(-1.00/(Y('+str(indexH+1)+')+1.89*Y('+str(indexH2+1)+')))\n'

                    # Photoreactions due to X-ray induced secondary photons (Lyman-alpha from excited H)
                    if reactants[i].count('XRLYA') == 1:
                        matrixString += '*Y('+str(indexH+1)+')*ZETA_H'
#                        additionalString += TO BE ADDED!!!
#                        additionalString += TO BE ADDED!!!

                    # Photoreactions due to X-ray induced secondary photons (Lyman-Werner from excited H2)
                    if reactants[i].count('XRPHOT') == 1:
                        matrixString += '*Y('+str(indexH2+1)+')*ZETA_H2'
#                        additionalString += TO BE ADDED!!!
#                        additionalString += TO BE ADDED!!!

            if matrixString != '':
                matrixString = '      PD('+str(formatCode % (n+1))+','+str(formatCode % (m+1))+') = '+matrixString+'\n'
                output.write(matrixString)

    # If X-ray reactions are present, write their additional partial derivative terms
    if xrayReactions:
        output.write(additionalString)

    # Write the function footer
    fileFooter = '\n      RETURN\n      END\n'
    output.write(fileFooter)
    output.close()

# Write the Jacobian matrix file in F77 language format
def write_jac_f77(fileName, speciesList, reactants, products):
    nSpecies = len(speciesList)
    nReactions = len(reactants)
    output = open(fileName, mode='w')

    # Determine if X-ray reactions are present in the chemical network
    if sum([reactantList.count('XRAY')+reactantList.count('XRSEC') for reactantList in reactants]) > 0:
        xrayReactions = True
    else:
        xrayReactions = False

    # Find the index numbers for H, H2 and He in the species list
    if 'H'  in speciesList: indexH  = speciesList.index('H')
    if 'H2' in speciesList: indexH2 = speciesList.index('H2')
    if 'He' in speciesList: indexHe = speciesList.index('He')

    # Specify the appropriate format code to represent the matrix indices
    if nSpecies >= 1000:
        formatCode = '%4i'
    elif nSpecies >=100:
        formatCode = '%3i'
    elif nSpecies >= 10:
        formatCode = '%2i'
    else:
        formatCode = '%i'

    # Write the comments and function header
    if xrayReactions:
        fileHeader = "      SUBROUTINE JAC(NEQ,T,Y,ML,MU,PD,NROWPD)\n      INCLUDE 'header.f'\n      INTEGER          NEQ,ML,MU,NROWPD\n      DOUBLE PRECISION T,D,Y,PD\n      DOUBLE PRECISION ZETA_H,ZETA_H2,ZETA_HE\n      DIMENSION        NEQ(*),Y(*),PD(NROWPD,*)\n\nC     Set D to the gas density for use in the matrix\n      D=DENS(DEPTH)\n"
    else:
        fileHeader = "      SUBROUTINE JAC(NEQ,T,Y,ML,MU,PD,NROWPD)\n      INCLUDE 'header.f'\n      INTEGER          NEQ,ML,MU,NROWPD\n      DOUBLE PRECISION T,D,Y,PD\n      DIMENSION        NEQ(*),Y(*),PD(NROWPD,*)\n\nC     Set D to the gas density for use in the matrix\n      D=DENS(DEPTH)\n"
    output.write(fileHeader)

    # If X-ray reactions are present, write the additional terms needed to calculate their rate partial derivatives
    if xrayReactions:
        output.write(xray_parameters(speciesList, codeFormat='F77'))
        additionalString = '\nC     The additional partial derivative terms for the X-ray secondary ionization reactions begin here...\n'

    # Prepare and write the terms for each Jacobian matrix element
    output.write('\nC     The Jacobian matrix created by MakeRates begins here...\n')
    if useTtk:
        app.progressValue.set(0)
        Style().configure('Horizontal.TProgressbar', background='#dcdad5')
    for n in range(nSpecies):
        if useTtk:
            app.progressValue.set(100.0*float(n)/float(nSpecies))
            app.status.update_idletasks()
        species1 = speciesList[n]
        for m in range(nSpecies):
            species2 = speciesList[m]
            matrixString = ''

            # Loss terms for species1
            for i in range(nReactions):
                if reactants[i].count(species1) > 0 and reactants[i].count(species2) > 0:

                    # Gas-grain reaction rates include a density factor to account for the number density of grains
                    if is_gasgrain_reaction(reactants[i], products[i]):
                        matrixString += '-'+multiple(reactants[i].count(species1))+'K('+str(i+1)+')*D'
                        continue

                    matrixString += '-'+multiple(reactants[i].count(species1))+multiple(reactants[i].count(species1))+'K('+str(i+1)+')'

                    # Purely gas phase or grain surface reactions
                    for reactant in speciesList:
                        if reactant == species2:
                            for j in range(reactants[i].count(reactant)-1):
                                matrixString += '*Y('+str(speciesList.index(reactant)+1)+')*D'
                        else:
                            for j in range(reactants[i].count(reactant)):
                                matrixString += '*Y('+str(speciesList.index(reactant)+1)+')*D'

                    # Multiply by the electron density when electrons appear as a reactant
                    for j in range(reactants[i].count('e-')):
                        matrixString += '*X(1)*D'

                    # X-ray induced secondary ionization
                    if reactants[i].count('XRSEC') == 1:
                        if reactants[i].count('H') == 1:
                            matrixString = matrixString[:-len('-K('+str(i+1)+')')]
                            additionalString += '      PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH2+1))+') = PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH2+1))+')-K('+str(i+1)+')*Y('+str(indexH+1) +')*ZETA_H*(-1.89/(Y('+str(indexH+1)+')+1.89*Y('+str(indexH2+1)+')))\n'
                            additionalString += '      PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH+1)) +') = PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH+1)) +')-K('+str(i+1)+')*Y('+str(indexH2+1)+')*ZETA_H*(+1.89/(Y('+str(indexH+1)+')+1.89*Y('+str(indexH2+1)+')))\n'
                        elif reactants[i].count('H2') == 1:
                            matrixString = matrixString[:-len('-K('+str(i+1)+')')]
                            additionalString += '      PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH2+1))+') = PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH2+1))+')-K('+str(i+1)+')*Y('+str(indexH+1) +')*ZETA_H2*(+0.53/(Y('+str(indexH2+1)+')+0.53*Y('+str(indexH+1)+')))\n'
                            additionalString += '      PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH+1)) +') = PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH+1)) +')-K('+str(i+1)+')*Y('+str(indexH2+1)+')*ZETA_H2*(-0.53/(Y('+str(indexH2+1)+')+0.53*Y('+str(indexH+1)+')))\n'
                        elif reactants[i].count('He') == 1:
                            matrixString += '*ZETA_HE'
                        else:
                            matrixString += '*ZETA_H'
                            additionalString += '      PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH2+1))+') = PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH2+1))+')-K('+str(i+1)+')*Y('+str(n+1)+')*ZETA_H*(-1.89/(Y('+str(indexH+1)+')+1.89*Y('+str(indexH2+1)+')))\n'
                            additionalString += '      PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH+1)) +') = PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH+1)) +')-K('+str(i+1)+')*Y('+str(n+1)+')*ZETA_H*(-1.00/(Y('+str(indexH+1)+')+1.89*Y('+str(indexH2+1)+')))\n'

                    # Photoreactions due to X-ray induced secondary photons (Lyman-alpha from excited H)
                    if reactants[i].count('XRLYA') == 1:
                        matrixString += '*Y('+str(indexH+1)+')*ZETA_H'
#                        additionalString += TO BE ADDED!!!
#                        additionalString += TO BE ADDED!!!

                    # Photoreactions due to X-ray induced secondary photons (Lyman-Werner from excited H2)
                    if reactants[i].count('XRPHOT') == 1:
                        matrixString += '*Y('+str(indexH2+1)+')*ZETA_H2'
#                        additionalString += TO BE ADDED!!!
#                        additionalString += TO BE ADDED!!!

            # Formation terms for species1
            for i in range(nReactions):
                if products[i].count(species1) > 0 and reactants[i].count(species2) > 0:
                    matrixString += '+'+multiple(products[i].count(species1))+'K('+str(i+1)+')'

                    # Gas-grain reaction rates include a density factor to account for the number density of grains
                    if is_gasgrain_reaction(reactants[i], products[i]):
                        matrixString += '*D'
                        continue

                    # Purely gas phase or grain surface reactions
                    for reactant in speciesList:
                        if reactant == species2:
                            for j in range(reactants[i].count(reactant)-1):
                                matrixString += '*Y('+str(speciesList.index(reactant)+1)+')*D'
                        else:
                            for j in range(reactants[i].count(reactant)):
                                matrixString += '*Y('+str(speciesList.index(reactant)+1)+')*D'

                    # Multiply by the electron density when electrons appear as a reactant
                    for j in range(reactants[i].count('e-')):
                        matrixString += '*X(1)*D'

                    # X-ray induced secondary ionization
                    if reactants[i].count('XRSEC') == 1:
                        if reactants[i].count('H') == 1:
                            matrixString = matrixString[:-len('+K('+str(i+1)+')')]
                            additionalString += '      PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH2+1))+') = PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH2+1))+')+K('+str(i+1)+')*Y('+str(indexH+1) +')*ZETA_H*(-1.89/(Y('+str(indexH+1)+')+1.89*Y('+str(indexH2+1)+')))\n'
                            additionalString += '      PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH+1)) +') = PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH+1)) +')+K('+str(i+1)+')*Y('+str(indexH2+1)+')*ZETA_H*(+1.89/(Y('+str(indexH+1)+')+1.89*Y('+str(indexH2+1)+')))\n'
                        elif reactants[i].count('H2') == 1:
                            matrixString = matrixString[:-len('+K['+str(i)+']')]
                            additionalString += '      PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH2+1))+') = PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH2+1))+')+K('+str(i+1)+')*Y('+str(indexH+1) +')*ZETA_H2*(+0.53/(Y('+str(indexH2+1)+')+0.53*Y('+str(indexH+1)+')))\n'
                            additionalString += '      PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH+1)) +') = PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH+1)) +')+K('+str(i+1)+')*Y('+str(indexH2+1)+')*ZETA_H2*(-0.53/(Y('+str(indexH2+1)+')+0.53*Y('+str(indexH+1)+')))\n'
                        elif reactants[i].count('He') == 1:
                            matrixString += '*ZETA_HE'
                        else:
                            matrixString += '*ZETA_H'
                            additionalString += '      PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH2+1))+') = PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH2+1))+')+K('+str(i+1)+')*Y('+str(m+1)+']*ZETA_H*(-1.89/(Y('+str(indexH+1)+')+1.89*Y('+str(indexH2+1)+')))\n'
                            additionalString += '      PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH+1)) +') = PD('+str(formatCode % (n+1))+','+str(formatCode % (indexH+1)) +')+K('+str(i+1)+')*Y('+str(m+1)+']*ZETA_H*(-1.00/(Y('+str(indexH+1)+')+1.89*Y('+str(indexH2+1)+')))\n'

                    # Photoreactions due to X-ray induced secondary photons (Lyman-alpha from excited H)
                    if reactants[i].count('XRLYA') == 1:
                        matrixString += '*Y('+str(indexH+1)+')*ZETA_H'
#                        additionalString += TO BE ADDED!!!
#                        additionalString += TO BE ADDED!!!

                    # Photoreactions due to X-ray induced secondary photons (Lyman-Werner from excited H2)
                    if reactants[i].count('XRPHOT') == 1:
                        matrixString += '*Y('+str(indexH2+1)+')*ZETA_H2'
#                        additionalString += TO BE ADDED!!!
#                        additionalString += TO BE ADDED!!!

            if matrixString != '':
                matrixString = '      PD('+str(formatCode % (n+1))+','+str(formatCode % (m+1))+') = '+matrixString+'\n'
                matrixString = truncate_line(matrixString, continuationCode='     &              ')
                output.write(matrixString)

    # If X-ray reactions are present, write their additional partial derivative terms
    if xrayReactions:
        additionalString = truncate_line(additionalString, continuationCode='     &              ')
        output.write(additionalString)

    # Write the function footer
    fileFooter = '\n      RETURN\n      END\n'
    output.write(fileFooter)
    output.close()


# --- Application declaration begins here --- #


class Application(Frame):
    def __init__(self, master=None):
        Frame.__init__(self, master)
        if useTtk:
            Style().theme_use('clam')
        self.grid()
        self.createWidgets()
        self.setupWidgets()

    def reactionFilePrompt(self):
        filename = tkFileDialog.askopenfilename(defaultextension='.dat', title='Choose reaction file')
        if filename:
            self.reactionFile.set(filename)

    def speciesFilePrompt(self):
        filename = tkFileDialog.askopenfilename(defaultextension='.dat', title='Choose species file')
        if filename:
            self.speciesFile.set(filename)

    def statusMessage(self, messageString, replace=False, error=False):
        if replace:
            self.status.delete('1.0', 'end')
        if error:
            self.status.insert('end', messageString, ('Error'))
            if useTtk:
                Style().configure('Horizontal.TProgressbar', background='red')
                self.progressValue.set(100)
        else:
            self.status.insert('end', messageString, ('Normal'))
        self.status.see('end')
        self.status.update_idletasks()

    def executeCode(self):

        # Check that the specified input files exist and that the supplied options are valid
        outputPrefix = self.outputPrefix.get()
        if outputPrefix == '<None>': outputPrefix = ''

        reactionFile = self.reactionFile.get()
        if reactionFile == '<None>' or reactionFile == '':
            self.statusMessage('ERROR! A input reaction file must be specified.', replace=True, error=True)
            return
        if not os.path.isfile(reactionFile):
            self.statusMessage('ERROR! Specified reaction file ('+reactionFile+') does not exist.', replace=True, error=True)
            return
        if reactionFile == outputPrefix+'rates.dat':
            self.outputPrefix.set('output-')
            outputPrefix = 'output-'

        speciesFile = self.speciesFile.get()
        if speciesFile == '<None>': speciesFile = ''
        if speciesFile != '':
            if not os.path.isfile(speciesFile):
                self.statusMessage('ERROR! Specified species file ('+speciesFile+') does not exist.', replace=True, error=True)
                return
            if speciesFile == outputPrefix+'species.dat':
                self.outputPrefix.set('output-')
                outputPrefix = 'output-'

        sortSpecies = self.sortSpecies.get()
        if sortSpecies == 1: sortSpecies = True
        else: sortSpecies = False

        logForm = self.logFormat.get()
        if logForm == 1: logForm = True
        else: logForm = False

        fileFormat = self.fileFormatList.curselection()
        if len(fileFormat) > 0: fileFormat = fileFormat[0]
        if   fileFormat == '0': fileFormat = 'rate05'
        elif fileFormat == '1': fileFormat = 'rate99'
        elif fileFormat == '2': fileFormat = 'rate95'

        codeFormat = self.codeFormatList.curselection()
        if len(codeFormat) > 0: codeFormat = codeFormat[0]
        if   codeFormat == '0': codeFormat = 'C'
        elif codeFormat == '1': codeFormat = 'F90'
        elif codeFormat == '2': codeFormat = 'F77'

        # Read the reactants, products, Arrhenius equation parameters and measurement labels for each reaction
        self.statusMessage('Reading reaction file...', replace=True)
        nReactions, reactants, products, alpha, beta, gamma, labels = read_reaction_file(reactionFile)
        self.statusMessage('\n')

        # Read the name, abundance and molecular mass for each species
        if speciesFile != '':
            self.statusMessage('Reading species file...')
            nSpecies, speciesList, abundanceList, massList = read_species_file(speciesFile)
            self.statusMessage('\n')

            # Find the total number and full list of reactions containing only these species
            self.statusMessage('\nFinding all reactions involving these species...')
            nReactions, reactants, products, alpha, beta, gamma, labels = find_all_reactions(speciesList, reactants, products, alpha, beta, gamma, labels)
            self.statusMessage('\n')

        # Find the total number and full list of unique species contained in the reactions
        else:
            self.statusMessage('\nFinding all species involved in these reactions...')
            nSpecies, speciesList = find_all_species(reactants, products)
            abundanceList = [float(0) for i in range(nSpecies)]
            self.statusMessage('\n')

        # Discard species containing more than a certain number of atoms (not counting H or D)
        speciesSizeLimit = 0
        if speciesSizeLimit > 0:
            if speciesSizeLimit == 1:
                self.statusMessage('\nDiscarding species containing more than '+str(speciesSizeLimit)+' atom \n(not counting H or D atoms)...')
            else:
                self.statusMessage('\nDiscarding species containing more than '+str(speciesSizeLimit)+' atoms\n(not counting H or D atoms)...')

            massList, constituentList, elementList = find_constituents(speciesList)
            index = elementList.index('D') + 1
            for i in range(nSpecies-1,-1,-1):
                if sum(constituentList[i][index:]) > speciesSizeLimit:
                    del(speciesList[i], massList[i], constituentList[i])

            # Update the total number of species
            nSpecies = len(speciesList)

            # Find the total number and full list of reactions containing only this new reduced set of species
            nReactions, reactants, products, alpha, beta, gamma, labels = find_all_reactions(speciesList, reactants, products, alpha, beta, gamma, labels)
            self.statusMessage('\n')

        self.statusMessage('\nNumber of reactions: '+str(nReactions)+'\nNumber of species: '+str(nSpecies)+'\n')

        # Check for "orphan" species that are either never formed or never destroyed
        self.statusMessage('\nChecking for species without formation/destruction reactions...')
        nFormation, nDestruction, missingList = check_orphan_species(speciesList, reactants, products)
        self.statusMessage('\n')

        # Sort the species first by number of destruction reactions, then by number of formation reactions
        if sortSpecies:
            self.statusMessage('\nSorting the species by number of formation reactions...')
            speciesList, abundanceList = sort_species(speciesList, abundanceList, nFormation, nDestruction)
            self.statusMessage('\n')

        # Calculate the molecular mass and elemental constituents of each species
        self.statusMessage('\nCalculating molecular masses and elemental constituents...')
        massList, constituentList, elementList = find_constituents(speciesList)
        self.statusMessage('\n')

        # Check that all elemental constituents and the net charge are conserved in each reaction
        self.statusMessage('\nChecking element and charge conservation for each reaction...')
        failed = check_conservation(reactants, products, speciesList, elementList, constituentList)
        if failed: return
        self.statusMessage('\n')

        # Create the species file
        self.statusMessage('\nWriting species file in '+fileFormat.title()+' format...\n')
        filename = outputPrefix+'species.dat'
        write_species(filename, speciesList, abundanceList, massList, fileFormat=fileFormat)

        # Create the reaction file
        self.statusMessage('Writing reaction file in '+fileFormat.title()+' format...\n')
        filename = outputPrefix+'rates.dat'
        write_reactions(filename, reactants, products, alpha, beta, gamma, labels, fileFormat=fileFormat)

        # Write the ODEs in the appropriate language format
        if codeFormat == 'C':
            self.statusMessage('Writing system of ODEs in C format...\n')
            filename = outputPrefix+'odes.c'
            write_odes_c(filename, speciesList, constituentList, reactants, products, logForm=logForm)
        if codeFormat == 'F90':
            self.statusMessage('Writing system of ODEs in F90 format...\n')
            filename = outputPrefix+'odes.f90'
            write_odes_f90(filename, speciesList, constituentList, reactants, products)
        if codeFormat == 'F77':
            self.statusMessage('Writing system of ODEs in F77 format...\n')
            filename = outputPrefix+'odes.f'
            write_odes_f77(filename, speciesList, constituentList, reactants, products)

        # Write the Jacobian matrix in the appropriate language format
        if codeFormat == 'C':
            self.statusMessage('Writing Jacobian matrix in C format...\n')
            filename = outputPrefix+'jacobian.c'
            write_jac_c(filename, speciesList, reactants, products, logForm=logForm)
        if codeFormat == 'F90':
            self.statusMessage('Writing Jacobian matrix in F90 format...\n')
            filename = outputPrefix+'jac.f90'
            write_jac_f90(filename, speciesList, reactants, products)
        if codeFormat == 'F77':
            self.statusMessage('Writing Jacobian matrix in F77 format...\n')
            filename = outputPrefix+'jac.f'
            write_jac_f77(filename, speciesList, reactants, products)

        if useTtk:
            Style().configure('Horizontal.TProgressbar', background='#00cc00')
            self.progressValue.set(100)
        self.statusMessage('\nFinished!')
        return

    def createWidgets(self):

        # Create a frame to contain the filename entry fields
        self.filenameFrame = Frame(self)
        self.filenameFrame.grid(row=0, column=0, columnspan=2, padx=5, pady=5)

        # Create the reaction filename entry field and selection dialog button
        self.reactionFileLabel = Label(self.filenameFrame, text='Reaction File: ')
        self.reactionFileLabel.grid(row=0, column=0, sticky='E')
        self.reactionFileEntry = Entry(self.filenameFrame)
        self.reactionFileEntry.grid(row=0, column=1)
        if useTtk:
            Style().configure('BrowseFile.TButton', font='helvetica 11', height=1, width=2, padding=-4)
            self.reactionFileSelect = Button(self.filenameFrame, command=self.reactionFilePrompt, text='+',style='BrowseFile.TButton')
        else:
            self.reactionFileSelect = Button(self.filenameFrame, command=self.reactionFilePrompt, bitmap='gray50', height=6, width=6, foreground='#0088ff', background='#0088ff', activeforeground='#22aaff', activebackground='#22aaff')
        self.reactionFileSelect.grid(row=0, column=2, padx=1)

        # Create the species filename entry field and selection dialog button
        self.speciesFileLabel = Label(self.filenameFrame, text='Species File: ')
        self.speciesFileLabel.grid(row=1, column=0, sticky='E')
        self.speciesFileEntry = Entry(self.filenameFrame)
        self.speciesFileEntry.grid(row=1, column=1)
        if useTtk:
            self.speciesFileSelect = Button(self.filenameFrame, command=self.speciesFilePrompt, text='+', style='BrowseFile.TButton')
        else:
            self.speciesFileSelect = Button(self.filenameFrame, command=self.speciesFilePrompt, bitmap='gray50', height=6, width=6, foreground='#0088ff', background='#0088ff', activeforeground='#22aaff', activebackground='#22aaff')
        self.speciesFileSelect.grid(row=1, column=2, padx=1)

        # Create the entry field for the output file prefix
        self.outputPrefixLabel = Label(self.filenameFrame, text='Output Prefix: ')
        self.outputPrefixLabel.grid(row=2, column=0, sticky='E')
        self.outputPrefixEntry = Entry(self.filenameFrame)
        self.outputPrefixEntry.grid(row=2, column=1)

        # Create the checkbox to enable sorting of the species (by number of formation reactions)
        self.sortButton = Checkbutton(self, text='Sort Species')
        self.sortButton.grid(row=3, column=0)

        # Create the checkbox to enable writing of the ODEs in logarithmic form
        self.logButton = Checkbutton(self, text='Log Format')
        self.logButton.grid(row=3, column=1)

        # Create the list of output file formats (Rate05, Rate99, Rate95)
        self.fileFormatList = Listbox(self, exportselection=0, height=3, width=8)
        self.fileFormatList.grid(row=4, column=0)

        # Create the list of output code formats (C, F90, F77)
        self.codeFormatList = Listbox(self, exportselection=0, height=3, width=8)
        self.codeFormatList.grid(row=4, column=1)

        # Create the button to execute the code
        if useTtk:
            Style().configure('Go.TButton', foreground='#00cc00', activeforeground='#00cc00', width=7, font='helvetica 9 bold')
            self.go = Button(self, text='Go!', command=self.executeCode, style='Go.TButton')
        else:
            self.go = Button(self, text='Go!', command=self.executeCode, foreground='#00cc00', activeforeground='#00cc00')
        self.go.grid(row=5, column=0, pady=5)

        # Create the button to quit the application
        if useTtk:
            Style().configure('Quit.TButton', foreground='red', activeforeground='red', width=7, font='helvetica 9 bold')
            self.quit = Button(self, text='Exit', command=self.quit, style='Quit.TButton')
        else:
            self.quit = Button(self, text='Exit', command=self.quit, foreground='red', activeforeground='red')
        self.quit.grid(row=5, column=1, pady=5)

        # Create a frame to contain the status message area and associated scrollbar
        self.messageFrame = Frame(self)
        self.messageFrame.grid(row=6, column=0, columnspan=2, padx=5, pady=5)

        # Create the status message area
        self.status = Text(self.messageFrame, exportselection=0, height=10, width=40, undo=False, wrap=WORD, font=tkFont.Font(family='Helvetica', size=10))
        self.status.grid(row=0, column=0, ipadx=5, ipady=5, sticky=N+E+S+W)

        # Create a vertical scrollbar and connect it to the status message text area
        self.scrollbar = Scrollbar(self.messageFrame, orient=VERTICAL, command=self.status.yview)
        self.scrollbar.grid(row=0, column=1, sticky=N+S)
        self.status['yscrollcommand'] = self.scrollbar.set

        # Create a horizontal progress bar to indicate the status of each task
        if useTtk:
            Style().configure('Horizontal.TProgressbar')
            self.progressbar = Progressbar(self.messageFrame, orient=HORIZONTAL, mode='determinate', style='Horizontal.TProgressbar')
            self.progressbar.grid(row=1, column=0, padx=1, sticky=N+E+S+W)

    def setupWidgets(self):

        # Specify the reaction file name variable
        self.reactionFile = StringVar()
        # Set it to the default value
        self.reactionFile.set(reactionFile)
        # Tell the associated entry widget to watch this variable
        self.reactionFileEntry['textvariable'] = self.reactionFile

        # Specify the species file name variable
        self.speciesFile = StringVar()
        # Set it to the default value
        self.speciesFile.set(speciesFile)
        # Tell the associated entry widget to watch this variable
        self.speciesFileEntry['textvariable'] = self.speciesFile

        # Specify the output file prefix variable
        self.outputPrefix = StringVar()
        # Set it to the default value
        self.outputPrefix.set(outputPrefix)
        # Tell the associated entry widget to watch this variable
        self.outputPrefixEntry['textvariable'] = self.outputPrefix

        # Specify the sort species flag variable
        self.sortSpecies = IntVar()
        # Set it to the default value
        self.sortSpecies.set(sortSpecies)
        # Tell the associated checkbutton widget to watch this variable
        self.sortButton['variable'] = self.sortSpecies

        # Specify the logarithmic format flag variable
        self.logFormat = IntVar()
        # Set it to the default value
        self.logFormat.set(logForm)
        # Tell the associated checkbutton widget to watch this variable
        self.logButton['variable'] = self.logFormat

        # Specify the file format options variable
        self.fileFormat = StringVar()
        # Specify the available values
        self.fileFormat.set('Rate05 Rate99 Rate95')
        # Tell the associated listbox widget to watch this variable
        self.fileFormatList['listvariable'] = self.fileFormat
        # Select the default entry in the list
        self.fileFormatList.selection_set(['Rate05','Rate99','Rate95'].index(fileFormat))

        # Specify the code format options variable
        self.codeFormat = StringVar()
        # Specify the available values
        self.codeFormat.set('C F90 F77')
        # Tell the associated listbox widget to watch this variable
        self.codeFormatList['listvariable'] = self.codeFormat
        # Select the default entry in the list
        self.codeFormatList.selection_set(['C','F90','F77'].index(codeFormat))

        if useTtk:
            # Specify the progress bar status variable
            self.progressValue = IntVar()
            # Set it to the default value
            self.progressValue.set(0)
            # Tell the associated progressbar widget to watch this variable
            self.progressbar['variable'] = self.progressValue

        # Specify the text style tags for status messages
        self.status.tag_config('Error', foreground='red')
        self.status.tag_config('Normal', foreground='black')
        # Specify the welcome message to display at start
        self.statusMessage('Welcome to MakeRates.')


# --- Global code begins here --- #


# Read and verify the command line keyword options (if any have been specified)
usageString = 'MakeRates.py [reactionFile=ReactionFileName] [speciesFile=SpeciesFileName] [outputPrefix=OutputFilePrefix] [sortSpecies=True|False] [logForm=True|False] [fileFormat=Rate95|Rate99|Rate05] [codeFormat=F77|F90|C]\n'
keywordNames = ['reactionFile=','speciesFile=','outputPrefix=','sortSpecies=','logForm=','fileFormat=','codeFormat=']
keywordValue = [reactionFile,speciesFile,outputPrefix,('True'if sortSpecies else 'False'),('True'if logForm else 'False'),fileFormat,codeFormat]

if len(sys.argv) > 1 and sys.argv[1].lower() == '--help':
    sys.exit('Usage:\n\n'+usageString)

if len(sys.argv) > 1 and sys.argv[1].lower() == '--version':
    sys.exit('MakeRates version '+versionString)

for n in range(1,len(sys.argv)):
    if sum([sys.argv[n].lower().count(keyword.lower()) for keyword in keywordNames]) == 0:
        sys.exit('\nERROR! Unrecognised keyword: '+sys.argv[n]+'\n\nUsage:\n'+usageString)

    for i, keyword in enumerate(keywordNames):
        if sys.argv[n].lower().count(keyword.lower()) != 0:
            index = sys.argv[n].lower().index(keyword.lower())+len(keyword)
            keywordValue[i] = sys.argv[n][index:].strip()

if keywordValue[0] != '':
    reactionFile = keywordValue[0]
else:
    reactionFile = '<None>'

if keywordValue[1] != '':
    speciesFile = keywordValue[1]
else:
    speciesFile = '<None>'

if keywordValue[2] != '':
    outputPrefix = keywordValue[2]
else:
    outputPrefix = '<None>'

if keywordValue[3].title() == 'True' or keywordValue[3].title() == 'False':
    sortSpecies = keywordValue[3].title() == 'True'
else:
    sys.exit('\nERROR! Unrecognised option for keyword '+keywordNames[3]+': '+keywordValue[3]+'\n\nUsage: '+usageString)

if keywordValue[4].title() == 'True' or keywordValue[4].title() == 'False':
    logForm = keywordValue[4].title() == 'True'
else:
    sys.exit('\nERROR! Unrecognised option for keyword '+keywordNames[4]+': '+keywordValue[4]+'\n\nUsage: '+usageString)

if keywordValue[5].title() == 'Rate95' or keywordValue[5].title() == 'Rate99' or keywordValue[5].title() == 'Rate05':
    fileFormat = keywordValue[5].title()
else:
    sys.exit('\nERROR! Unrecognised option for keyword '+keywordNames[5]+': '+keywordValue[5]+'\n\nUsage: '+usageString)

if keywordValue[6].upper() == 'F77' or keywordValue[6].upper() == 'F90' or keywordValue[6].upper() == 'C':
    codeFormat = keywordValue[6].upper()
else:
    sys.exit('\nERROR! Unrecognised option for keyword '+keywordNames[6]+': '+keywordValue[6]+'\n\nUsage: '+usageString)

# Start the Tkinter GUI application, if supported
if useTkinter:
    root = Tk()
    app = Application(master=root)
    app.master.title('MakeRates')
    app.mainloop()
    root.destroy()

# Otherwise, use command line print statements to inform the user of progress
else:
    if reactionFile != '':
        if not os.path.isfile(reactionFile):
            sys.exit('\nERROR! Specified reaction file '+reactionFile+' does not exist\n')
    else:
        sys.exit('\nERROR! An input reaction file must be specified\n\nUsage: '+usageString)

    if speciesFile != '':
        if not os.path.isfile(speciesFile):
            sys.exit('\nERROR! Specified species file '+speciesFile+' does not exist\n')

    # Read the reactants, products, Arrhenius equation parameters and measurement labels for each reaction
    print '\nReading reaction file...'
    nReactions, reactants, products, alpha, beta, gamma, labels = read_reaction_file(reactionFile)

    # Read the name, abundance and molecular mass for each species
    if speciesFile != '':
        print 'Reading species file...'
        nSpecies, speciesList, abundanceList, massList = read_species_file(speciesFile)

    # Find the total number and full list of reactions containing only these species
        print '\nFinding all reactions involving these species...'
        nReactions, reactants, products, alpha, beta, gamma, labels = find_all_reactions(speciesList, reactants, products, alpha, beta, gamma, labels)

    # Find the total number and full list of unique species contained in the reactions
    else:
        print '\nFinding all species involved in the reactions...'
        nSpecies, speciesList = find_all_species(reactants, products)
        abundanceList = [float(0) for i in range(nSpecies)]

    print '\nNumber of reactions:',nReactions
    print 'Number of species:',nSpecies

    # Check for "orphan" species that are either never formed or never destroyed
    print '\nChecking for species without formation/destruction reactions...'
    nFormation, nDestruction, missingList = check_orphan_species(speciesList, reactants, products)

    # Sort the species first by number of destruction reactions, then by number of formation reactions
    if sortSpecies:
        print '\nSorting the species by number of formation reactions...'
        speciesList, abundanceList = sort_species(speciesList, abundanceList, nFormation, nDestruction)

    # Calculate the molecular mass and elemental constituents of each species
    print '\nCalculating molecular masses and elemental constituents...'
    massList, constituentList, elementList = find_constituents(speciesList)

    # Check that all elemental constituents and the net charge are conserved in each reaction
    print '\nChecking element and charge conservation for each reaction...'
    status = check_conservation(reactants, products, speciesList, elementList, constituentList)

    # Create the species file
    print '\nWriting species file in '+fileFormat.title()+' format...'
    filename = outputPrefix+'species.dat'
    write_species(filename, speciesList, abundanceList, massList, fileFormat=fileFormat)

    # Create the reaction file
    print 'Writing reaction file in '+fileFormat.title()+' format...'
    filename = outputPrefix+'rates.dat'
    write_reactions(filename, reactants, products, alpha, beta, gamma, labels, fileFormat=fileFormat)

    # Write the ODEs in the appropriate language format
    if codeFormat == 'C':
        print 'Writing system of ODEs in C format...'
        filename = outputPrefix+'odes.c'
        write_odes_c(filename, speciesList, constituentList, reactants, products, logForm=logForm)
    if codeFormat == 'F90':
        print 'Writing system of ODEs in F90 format...'
        filename = outputPrefix+'odes.f90'
        write_odes_f90(filename, speciesList, constituentList, reactants, products)
    if codeFormat == 'F77':
        print 'Writing system of ODEs in F77 format...'
        filename = outputPrefix+'odes.f'
        write_odes_f77(filename, speciesList, constituentList, reactants, products)

    # Write the Jacobian matrix in the appropriate language format
    if codeFormat == 'C':
        print 'Writing Jacobian matrix in C format...'
        filename = outputPrefix+'jacobian.c'
        write_jac_c(filename, speciesList, reactants, products, logScale=False)
    if codeFormat == 'F90':
        print 'Writing Jacobian matrix in F90 format...'
        filename = outputPrefix+'jac.f90'
        write_jac_f90(filename, speciesList, reactants, products)
    if codeFormat == 'F77':
        print 'Writing Jacobian matrix in F77 format...'
        filename = outputPrefix+'jac.f'
        write_jac_f77(filename, speciesList, reactants, products)

    print '\nFinished!'

# --- End of the code --- #
