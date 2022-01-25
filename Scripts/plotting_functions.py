###############################################################################
# Functions to read and format uclpdr outputs
# Based on script by Tom A. Bell
# Amended by Jon Holdship
###############################################################################
import pandas as pd
import numpy as np
from glob import glob
from os import remove
#executes all the other read functions and obtains the values for each particle
#for a single ray
def create_output_dataframe(file_prefix,ray):
    '''inputs: file_prefix - path and first part of filename of models
    ray - integer for which ray to use for Av
    output: pandas dataframe of all cloud properties'''
    
    #create dataframe
    output_df=pd.DataFrame()

    #load cloud properties
    #nParticles, gasDensity, dustDensity, gasTemperature, dustTemperature, fuvFlux, xrayFlux=read_cloud_properties(file_prefix+"prop.out")
    property_names=["n_H","n_g","T_g","T_d","FUV","F_x"]
    nParticles,properties=read_cloud_properties(file_prefix+".prop.out")

    #get particle list and make column for it
    particles=np.arange(nParticles)
    particles=particles+1
    output_df["Particle"]=particles

    #then one column for each gas property
    for i,property_name in enumerate(property_names):
        output_df[property_name]=properties[i]

    #get particle avs for all rays
    nRays, particleCoordinates, visualExtinction = read_visual_extinctions(file_prefix+".av.out")
    #add just chosen ray to dataframe
    output_df["Av"]=[visualExtinction[i][ray] for i in range(nParticles)]

    #below is just repetitive code. That should be abstracted.
    #We load each of the uclpdr outputs and add columns to the dataframe for them
    #so we load abundances, loop through species adding one column per species abundance.

    #abundances
    nSpecies, speciesNames, fractionalAbundances = read_fractional_abundances(file_prefix+".abun.out")
    fractionalAbundances=np.asarray(fractionalAbundances)
    for i,speciesName in enumerate(speciesNames):
        output_df[speciesName+"_abun"]=fractionalAbundances[:,i]

    #cooling
    nCool, coolantNames, coolingRates=read_cooling_rates(file_prefix+".cool.out")
    coolingRates=np.asarray(coolingRates)
    for i, coolantName in enumerate(coolantNames):
        output_df[coolantName+"_cool"]=coolingRates[:,i]


    #heating
    nCool, coolantNames, coolingRates=read_cooling_rates(file_prefix+".heat.out")
    coolingRates=np.asarray(coolingRates)
    for i, coolantName in enumerate(coolantNames):
        output_df[coolantName+"_heat"]=coolingRates[:,i]    

    #opacities
    nLines, lineNames, lineOpacities=read_line_opacities(file_prefix+f".tau_{ray+1}.out")
    lineOpacities=np.asarray(lineOpacities)
    for i,lineName in enumerate(lineNames):
        output_df[lineName+"_tau"]=lineOpacities[:,i]

    #emissivities
    nLines, lineNames, localEmissivities = read_local_emissivities(file_prefix+".emis.out")
    localEmissivities=np.asarray(localEmissivities)
    for i,lineName in enumerate(lineNames):
        output_df[lineName+"_emis"]=localEmissivities[:,i]

    #population densities
    nLevels, populationNames, populationDensities = read_population_densities(file_prefix+".pop.out")
    populationDensities=np.asarray(populationDensities)
    for i,populationName in enumerate(populationNames):
        output_df[populationName+"_pop"]=populationDensities[:,i]

    #return complete dataframe
    return output_df

# Read the physical cloud properties
def read_cloud_properties(inputfile):
    inputfile=open(inputfile,"r")#.seek(0)
    data = inputfile.readline().split()
    gasDensity = [] ; dustDensity = [] ; gasTemperature = [] ; dustTemperature = [] ; fuvFlux = [] ; xrayFlux = []
    data = inputfile.readline().split()
    while data != []:
        gasDensity.append(number(data[1]))
        dustDensity.append(number(data[2]))
        gasTemperature.append(number(data[3]))
        dustTemperature.append(number(data[4]))
        fuvFlux.append(number(data[5]))
        xrayFlux.append(number(data[6]))
        data = inputfile.readline().split()
    nParticles = len(gasDensity)
    return nParticles, [gasDensity, dustDensity, gasTemperature, dustTemperature, fuvFlux, xrayFlux]

# Read the particle coordinates and the visual extinction along each ray to the PDR surface
def read_visual_extinctions(inputfile):
    inputfile=open(inputfile,"r")
    data = inputfile.readline().split()
    particleCoordinates = []; visualExtinction = []
    data = inputfile.readline().split()
    while data != []:
        nRays = len(data)-4
        particleCoordinates.append([number(data[i]) for i in range(1,4)])
        visualExtinction.append([number(data[i]) for i in range(4,len(data))])
        data = inputfile.readline().split()
    return nRays, particleCoordinates, visualExtinction

# Read the fractional abundances
def read_fractional_abundances(inputfile):
    inputfile=open(inputfile,"r")
    data = inputfile.readline().split()
    speciesName = [abundanceLabel[2:-1] for abundanceLabel in data[2:]]
    fractionalAbundance = []
    data = inputfile.readline().split()
    while data != []:
        nSpecies = len(data)-1
        fractionalAbundance.append([number(data[i]) for i in range(1,len(data))])
        data = inputfile.readline().split()
    return nSpecies, speciesName, fractionalAbundance

# Read the population densities
def read_population_densities(inputfile):
    inputfile=open(inputfile,"r")
    data = inputfile.readline().split()
    populationName = data[2:]
    populationDensity = []
    data = inputfile.readline().split()
    while data != []:
        nLevels = len(data)-1
        populationDensity.append([number(data[i]) for i in range(1,len(data))])
        data = inputfile.readline().split()
    return nLevels, populationName, populationDensity

# Read the local emissivities
def read_local_emissivities(inputfile):
    inputfile=open(inputfile,"r")
    data = inputfile.readline()[11:]
    lineName = []
    while data != '':
        lineName.append(r''+data[:13].strip())
        data = data[13:]
    for i in range(len(lineName)):
        if '[HI]' in lineName[i]:
            lineName[i] = r'Ly$\alpha$'
            continue
        lineName[i] = lineName[i].replace('[','$[$')
        lineName[i] = lineName[i].replace(']','$]$')
        lineName[i] = lineName[i].replace('III','{\sc\,iii}')
        lineName[i] = lineName[i].replace('II','{\sc\,ii}')
        lineName[i] = lineName[i].replace('I','{\sc\,i}')
        #lineName[i] = lineName[i].replace('um','$\upmu$m')
        lineName[i] = lineName[i].replace('A','\AA')
        lineName[i] = lineName[i].replace('-','--')
    localEmissivity = []
    data = inputfile.readline().split()
    while data != []:
        nLines = len(data)-1
        localEmissivity.append([number(data[i]) for i in range(1,len(data))])
        data = inputfile.readline().split()
    return nLines, lineName, localEmissivity

# Read the line opacities
def read_line_opacities(inputfile):
    inputfile=open(inputfile,"r")
    data = inputfile.readline()[11:]
    lineName = []
    while data != '':
        lineName.append(r''+data[:13].strip())
        data = data[13:]
    for i in range(len(lineName)):
        if '[HI]' in lineName[i]:
            lineName[i] = r'Ly$\alpha$'
            continue
        lineName[i] = lineName[i].replace('[','$[$')
        lineName[i] = lineName[i].replace(']','$]$')
       # lineName[i] = lineName[i].replace('III','{\sc\,iii}')
      #  lineName[i] = lineName[i].replace('II','{\sc\,ii}')
       # lineName[i] = lineName[i].replace('I','{\sc\,i}')
        #lineName[i] = lineName[i].replace('um','$\upmu$m')
        lineName[i] = lineName[i].replace('A','\AA')
        lineName[i] = lineName[i].replace('-','--')
    lineOpacity = []
    data = inputfile.readline().split()
    while data != []:
        nLines = len(data)-1
        lineOpacity.append([number(data[i]) for i in range(1,len(data))])
        data = inputfile.readline().split()
    return nLines, lineName, lineOpacity

# Read the cooling rates
def read_cooling_rates(inputfile):
    inputfile=open(inputfile,"r")
    data = inputfile.readline()[11:]
    coolantName = []
    while data != '':
        coolantName.append(r''+data[:13].strip())
        data = data[13:]
    for i in range(len(coolantName)):
        if '[HI]' in coolantName[i]:
            coolantName[i] = r'Ly$\alpha$'
            continue
        coolantName[i] = convert_species_name(coolantName[i])
        coolantName[i] = coolantName[i].replace('[','')
        coolantName[i] = coolantName[i].replace(']','')
        #coolantName[i] = coolantName[i].replace('III','{\sc\,iii}')
       # coolantName[i] = coolantName[i].replace('II','$\\sc{ii}$')
       # coolantName[i] = coolantName[i].replace('I','{\sc\,i}')
    coolingRate = []
    data = inputfile.readline().split()
    while data != []:
        nCool = len(data)-1
        coolingRate.append([number(data[i]) for i in range(1,len(data))])
        data = inputfile.readline().split()
    coolantName[-1] = r'Total'
    return nCool, coolantName, coolingRate

# Read the heating rates
def read_heating_rates(inputfile):
    inputfile=open(inputfile,"r")
    data = inputfile.readline().split()
    shortHeatingName = [r'Grain P.E.',r'H$_2$ Form',r'H$_2^*$ Pump',r'H$_2$ P.D.',r'Carbon P.I.',r'Cosmic-Rays',r'Turbulence',r'Chemistry',r'Gas-Grain',r'Coulomb',r'Total']
    heatingName = [r'Photoelectric',r'H$_2$ Formation',r'H$_2^*$ Pumping',r'H$_2$ Photodiss',r'C Ionization',r'Cosmic-Rays',r'Turbulence',r'Chemistry',r'Gas-Grain',r'Coulomb',r'Total Rate']
    heatingRate = []
    data = inputfile.readline().split()
    while data != []:
        nHeat = len(data)-1
        heatingRate.append([number(data[i]) for i in range(1,len(data))])
        data = inputfile.readline().split()
    heatingName[-1] = r'Total Rate'
    return nHeat, heatingName, heatingRate

# Convert a string representation of a number to a float, handling the
# case where the exponent is larger than +/-99 and the 'E' is missing
def number(string):
    try:
        return float(string)
    except ValueError:
        if string.count('+') > 0: fixedString = string[:string[2:].index('+')+2] + 'E' + string[string[2:].index('+')+2:]
        if string.count('-') > 0: fixedString = string[:string[2:].index('-')+2] + 'E' + string[string[2:].index('-')+2:]
        return float(fixedString)

# Produce the appropriate LaTeX string for a species name
def convert_species_name(name):
    string = ''
    for i in range(len(name)):
        if name[i:i+2] == 'p-': string += r'p-' ; continue
        if name[i:i+2] == 'o-': string += r'o-' ; continue
        if name[i-1:i+1] == 'p-': continue
        if name[i-1:i+1] == 'o-': continue
        if name[i].isalpha(): string += name[i]
        if name[i].isdigit(): string += r'$_{'+name[i]+'}$'
        if name[i] == '+': string += r'$^{+}$'
        if name[i] == '-': string += r'$^{-}$'
    string.replace('$$','')
    return string

def self(x):
    return x

def abundance_rename(species_name):
    species_name=species_name.replace("_abun","")
    species_name=species_name.upper()
    species_name=species_name.replace("2","$^2$")
    species_name=species_name.replace("3","$^3$")
    species_name=species_name.replace("+","$^+$")
    return species_name

def emissivity_rename(name):
    name=name.replace("_emis","")
    name=name.replace("um","$\mu$m")
    return name

def clear_raw_files(file_prefix):
    files=glob(file_prefix+"*")
    print(files)
    for file in files:
        remove(file)