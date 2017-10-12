#! /usr/bin/python
#
import codecs
import math
import os
import time
import numpy
import matplotlib
import matplotlib.colors as colors
import matplotlib.pyplot as pyplot
import matplotlib.gridspec as gridspec
import matplotlib.animation as animation
from matplotlib.ticker import FormatStrFormatter, LogLocator, MultipleLocator

# Specify the path for the ASCII data files
path = ''

# Specify the prefix and suffix for the model output data files
prefix = 'rollig_V1'
suffix = '.%04d'

# Specify the range of iterations to include
iterationRange = [1,20]

# Specify which data to plot
fluxPlot        = False
densityPlot     = False
densTempPlot    = False
temperaturePlot = False
abundancePlot   = False
populationPlot  = False
emissivityPlot  = False
opacityPlot     = False
coolingPlot     = True
heatingPlot     = False
heatCoolPlot    = False
opRatioPlot     = False

# Specify whether to plot visual extinction (A_V, True) or total column density (N_H, False) for the x-axis
plotVersusAV = False
plotVersusNH = not plotVersusAV

# Specify the plot axis ranges and major/minor tickmark spacings
extinctionAxisRange  = [[1E-4,2E+1], [10.,1.0]]
columnAxisRange      = [[1E20,1E25], [10.,1.0]]
fluxAxisRange        = [[1E-5,1E+3], [10.,1.0]]
densityAxisRange     = [[1E+1,1E+6], [10.,1.0]]
temperatureAxisRange = [[1E+1,1E+5], [10.,1.0]]
abundanceAxisRange   = [[1E-7,2E+1], [10.,1.0]]
populationAxisRange  = [[1E-10,1E+0], [100,100]]
emissivityAxisRange  = [[1E-28,1E-18], [100,100]]
coolingAxisRange     = [[1E-28,1E-18], [100,100]]
heatingAxisRange     = [[1E-28,1E-18], [100,100]]
opacityAxisRange     = [[1E-3,1E+2], [10.,1.0]]
intensityAxisRange   = [[1E-8,1E-3], [10.,1.0]]
opDensityAxisRange   = [[1E-5,1E+5], [100,100]]
opRatioAxisRange     = [[+0.0,+6.0], [100,100]]

# Specify which ray ID to plot against
rayID = 0

# Specify the species whose abundances should be plotted
# speciesList = ['H','H2','C+','C','CO','e-'] # PDR
speciesList = ['H','H2','O','C','CO','e-','H+','O+','C+'] # XDR

# Specify the indices of the level populations to be plotted
# levelList = [0,1,2,3,4] # C+
levelList = [5,6,7,8,9] # O
# levelList = [10,11,12,13,14] # C
# levelList = [15,16,17,18,19,20] # CO

# Specify the indices of the emission lines to be plotted
lineList = [1,2,4,9,11] # Fine structure lines
# lineList = [1] # [CII]
# lineList = [2,3,4,5,6,7,8] # [OI]
# lineList = [9,10,11,12,13,14,15] # [CI]
# lineList = [16,17,18,19,20,21,24] # CO

# Specify the lines whose integrated intensities should be plotted
emissionList = [" 63$"," 146$"," 158$"," 370$"," 609$",] # Fine structure lines
# emissionList = ["(1--0)","(2--1)","(3--2)","(4--3)","(5--4)"] # CO rotational lines

# Specify if the figure is to be used in a presentation
# (i.e. with white labels and a transparent background)
useForPresentation = False

# No user changes needed beyond this point...

# Read the physical cloud properties
def read_cloud_properties(inputfile):
    inputfile.seek(0)
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
    return nParticles, gasDensity, dustDensity, gasTemperature, dustTemperature, fuvFlux, xrayFlux

# Read the particle coordinates and the visual extinction along each ray to the PDR surface
def read_visual_extinctions(inputfile):
    inputfile.seek(0)
    data = inputfile.readline().split()
    particleCoordinates = []; visualExtinction = []
    data = inputfile.readline().split()
    while data != []:
        nRays = len(data)-4
        particleCoordinates.append([number(data[i]) for i in range(1,4)])
        visualExtinction.append([number(data[i]) for i in range(4,len(data))])
        data = inputfile.readline().split()
    return nRays, particleCoordinates, visualExtinction

# Read the gas temperatures and thermal balance data for the current iteration
def read_temperatures(inputfile):
    inputfile.seek(0)
    data = inputfile.readline().split()
    oldTemperature = [] ; newTemperature = [] ; lowerLimit = [] ; upperLimit = [] ; heatingRate = [] ; coolingRate = [] ; rateDifference = []
    data = inputfile.readline().split()
    while data != []:
        oldTemperature.append(number(data[1]))
        newTemperature.append(number(data[2]))
        lowerLimit.append(number(data[3]))
        upperLimit.append(number(data[4]))
        heatingRate.append(number(data[5]))
        coolingRate.append(number(data[6]))
        rateDifference.append(number(data[7]))
        data = inputfile.readline().split()
    return oldTemperature, newTemperature, lowerLimit, upperLimit, heatingRate, coolingRate, rateDifference

# Read the fractional abundances
def read_fractional_abundances(inputfile):
    inputfile.seek(0)
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
    inputfile.seek(0)
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
    inputfile.seek(0)
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
        lineName[i] = lineName[i].replace('III','\,{\sc iii}')
        lineName[i] = lineName[i].replace('II','\,{\sc ii}')
        lineName[i] = lineName[i].replace('I','\,{\sc i}')
        lineName[i] = lineName[i].replace('um','$\upmu$m')
        lineName[i] = lineName[i].replace('A','\AA')
        lineName[i] = lineName[i].replace('-','--')
    localEmissivity = []
    data = inputfile.readline().split()
    while data != []:
        nLines = len(data)-1
        localEmissivity.append([number(data[i]) for i in range(1,len(data))])
        data = inputfile.readline().split()
    return nLines, lineName, localEmissivity

# Read the cooling rates
def read_cooling_rates(inputfile):
    inputfile.seek(0)
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
        coolantName[i] = coolantName[i].replace('III','\,{\sc iii}')
        coolantName[i] = coolantName[i].replace('II','\,{\sc ii}')
        coolantName[i] = coolantName[i].replace('I','\,{\sc i}')
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
    inputfile.seek(0)
    data = inputfile.readline().split()
#    shortHeatingName = [r'Grain P.E.',r'H$_2$ Form.',r'H$_2^*$ Pumping',r'H$_2$ Photodiss',r'C Ioniz.',r'Cosmic-Rays',r'Turbulence',r'Chemistry',r'Gas-Grain',r'Total']
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

# Create a plot of the desired values
def create_plot(xData, yData, xAxisRange, yAxisRange, xAxisLabel, yAxisLabel, dataLabel, plotPosition, subplotPositions, logScale=[False,False], tickLabels=[True,True], axisLabels=[True,True], addLegend=False, color=False, style=False, linewidth=1.0):
    subplot = pyplot.subplot(subplotPositions[plotPosition-1])

    # Create a list of colour and line styles for the datasets, if not specified in the call
    if not color:
        color = ['red','orange','green','blue','purple','red','orange','green','blue','purple','red','orange','green','blue','purple','red','orange','green','blue','purple','red','orange','green','blue','purple','red','orange','green','blue','purple','red','orange','green','blue','purple','red','orange','green','blue','purple','red','orange','green','blue','purple','red','orange','green','blue','purple','red','orange','green','blue','purple','red','orange','green','blue','purple']
    if not style:
        style = ['solid','solid','solid','solid','solid','dashed','dashed','dashed','dashed','dashed','dotted','dotted','dotted','dotted','dotted','dashdot','dashdot','dashdot','dashdot','dashdot','solid','solid','solid','solid','solid','dashed','dashed','dashed','dashed','dashed','dotted','dotted','dotted','dotted','dotted','dashdot','dashdot','dashdot','dashdot','dashdot','solid','solid','solid','solid','solid','dashed','dashed','dashed','dashed','dashed','dotted','dotted','dotted','dotted','dotted','dashdot','dashdot','dashdot','dashdot','dashdot']

    # Plot each dataset
    plottedLines = []
    for n in range(len(yData)):
        line, = subplot.plot(xData, yData[n], color=color[n], linestyle=style[n], linewidth=1.0, label=dataLabel[n])
        plottedLines.append(line)

    # Set the axes scales, ranges and tickmarks
    if logScale[0]: pyplot.xscale('log')
    if logScale[1]: pyplot.yscale('log')
    pyplot.axis([xAxisRange[0][0], xAxisRange[0][1], yAxisRange[0][0], yAxisRange[0][1]])
    if logScale[0]:
        subplot.xaxis.set_major_locator(LogLocator(base=xAxisRange[1][0]))
        subplot.xaxis.set_minor_locator(LogLocator(base=xAxisRange[1][0], subs=range(int(xAxisRange[1][0]/xAxisRange[1][1]))))
    else:
        subplot.xaxis.set_major_locator(MultipleLocator(xAxisRange[1][0]))
        subplot.xaxis.set_minor_locator(MultipleLocator(xAxisRange[1][1]))
    if logScale[1]:
        subplot.yaxis.set_major_locator(LogLocator(base=yAxisRange[1][0]))
        subplot.yaxis.set_minor_locator(LogLocator(base=yAxisRange[1][0], subs=range(int(yAxisRange[1][0]/yAxisRange[1][1]))))
    else:
        subplot.yaxis.set_major_locator(MultipleLocator(yAxisRange[1][0]))
        subplot.yaxis.set_minor_locator(MultipleLocator(yAxisRange[1][1]))

    # Add tick and axis labels, if requested
    if not tickLabels[0]:
        pyplot.tick_params(axis='x', labelbottom='off')
    if not tickLabels[1]:
        pyplot.tick_params(axis='y', labelleft='off')
    if axisLabels[0]:
        pyplot.xlabel(xAxisLabel)
    if axisLabels[1]:
        pyplot.ylabel(yAxisLabel)

    # Add a legend, is requested
    if addLegend:
        if len(yData) > 9:
            ncol = 2
            mode = 'expand'
        elif max([len(label) for label in dataLabel]) > 20:
            ncol = 2
            mode = 'expand'
        elif len(yData) > 2:
            ncol = 3
            mode = 'expand'
        else:
            ncol = len(yData)
            mode = None
        pyplot.legend(mode=mode, loc='best', ncol=ncol, frameon=False, prop={'size':10}, borderpad=0, labelspacing=0.1, handletextpad=0.1, borderaxespad=0.6, handlelength=2.4)

    return plottedLines, subplot

# Beginning of the main program
start = time.time()

# Read the various data files
print '\nReading cloud properties...'
dataFilename = prefix+'.prop.out'
input = codecs.open(path+dataFilename, encoding='utf-8', mode='r')
nParticles, gasDensity, dustDensity, gasTemperature, dustTemperature, fuvFlux, xrayFlux = read_cloud_properties(input)
input.close()

print 'Reading visual extinctions...'
dataFilename = prefix+'.av.out'
input = codecs.open(path+dataFilename, encoding='utf-8', mode='r')
nRays, particleCoordinates, visualExtinction = read_visual_extinctions(input)
input.close()

if temperaturePlot:
    temperatureArray = []
    print 'Reading temperature data...'
    dataFilename = prefix+'.temp.0001'
    input = codecs.open(path+dataFilename, encoding='utf-8', mode='r')
    oldTemperature, newTemperature, lowerLimit, upperLimit, heatingRate, coolingRate, rateDifference = read_temperatures(input)
    input.close()
    temperatureArray.append([oldTemperature, oldTemperature, [temperature*0.7 for temperature in oldTemperature], [temperature*1.3 for temperature in oldTemperature]])
    for iteration in range(iterationRange[0],iterationRange[1]+1):
        dataFilename = prefix+'.temp.%04d' % iteration
        input = codecs.open(path+dataFilename, encoding='utf-8', mode='r')
        oldTemperature, newTemperature, lowerLimit, upperLimit, heatingRate, coolingRate, rateDifference = read_temperatures(input)
        input.close()
        temperatureArray.append([oldTemperature, newTemperature, lowerLimit, upperLimit])

if abundancePlot:
    abundanceArray = []
    print 'Reading fractional abundances...'
    dataFilename = prefix+'.abun.ini'
    input = codecs.open(path+dataFilename, encoding='utf-8', mode='r')
    nSpecies, speciesName, fractionalAbundance = read_fractional_abundances(input)
    input.close()
    abundanceArray.append(fractionalAbundance)
    for iteration in range(iterationRange[0],iterationRange[1]+1):
        dataFilename = prefix+'.abun.%04d' % iteration
        input = codecs.open(path+dataFilename, encoding='utf-8', mode='r')
        nSpecies, speciesName, fractionalAbundance = read_fractional_abundances(input)
        input.close()
        abundanceArray.append(fractionalAbundance)

if populationPlot:
    populationArray = []
    print 'Reading population densities...'
    dataFilename = prefix+'.pop.ini'
    input = codecs.open(path+dataFilename, encoding='utf-8', mode='r')
    nLevels, populationName, populationDensity = read_population_densities(input)
    input.close()
    populationArray.append(populationDensity)
    for iteration in range(iterationRange[0],iterationRange[1]+1):
        dataFilename = prefix+'.pop.%04d' % iteration
        input = codecs.open(path+dataFilename, encoding='utf-8', mode='r')
        nLevels, populationName, populationDensity = read_population_densities(input)
        input.close()
        populationArray.append(populationDensity)

if emissivityPlot:
    emissivityArray = []
    print 'Reading line emissivities...'
    for iteration in range(iterationRange[0],iterationRange[1]+1):
        dataFilename = prefix+'.emis.%04d' % iteration
        input = codecs.open(path+dataFilename, encoding='utf-8', mode='r')
        nLines, lineName, localEmissivity = read_local_emissivities(input)
        input.close()
        emissivityArray.append(localEmissivity)

if coolingPlot:
    coolingArray = []
    print 'Reading cooling rates...'
    for iteration in range(iterationRange[0],iterationRange[1]+1):
        dataFilename = prefix+'.cool.%04d' % iteration
        input = codecs.open(path+dataFilename, encoding='utf-8', mode='r')
        nCool, coolantName, coolingRate = read_cooling_rates(input)
        input.close()
        coolingArray.append(coolingRate)

if heatingPlot or coolingPlot:
    heatingArray = []
    print 'Reading heating rates...'
    for iteration in range(iterationRange[0],iterationRange[1]+1):
        dataFilename = prefix+'.heat.%04d' % iteration
        input = codecs.open(path+dataFilename, encoding='utf-8', mode='r')
        nHeat, heatingName, heatingRate = read_heating_rates(input)
        input.close()
        heatingArray.append(heatingRate)

if heatCoolPlot:
    heatCoolArray = []
    print 'Reading heating/cooling rates...'
    for iteration in range(iterationRange[0],iterationRange[1]+1):
        dataFilename = prefix+'.temp.%04d' % iteration
        input = codecs.open(path+dataFilename, encoding='utf-8', mode='r')
        oldTemperature, newTemperature, lowerLimit, upperLimit, heatingRate, coolingRate, rateDifference = read_temperatures(input)
        input.close()
        heatCoolArray.append([heatingRate, coolingRate])

# Convert the visual extinction values (A_V) to total column densities (N_H)
totalColumn = [[A_V_ray/6.289E-22 for A_V_ray in A_V] for A_V in visualExtinction]

# Convert the FUV values from multiples of the Draine field to fluxes in erg cm^-2 s^-1
fuvFlux = [chi*2.74E-3 for chi in fuvFlux]

# Add gas-grain collisions to the list of coolants when it acts as a negative heating rate
if coolingPlot:
    index = heatingName.index('Gas-Grain')
    if min(heatingArray[:][:][index]) < 0:
        nCool += 1
        coolantName.insert(-1,r'Grain')
        for i in range(len(coolingArray)):
            for n in range(nParticles):
                coolingArray[i][n].insert(-1,-heatingArray[i][n][index])

# Specify the axis labels
avLabel   = r'$A_V$ (mag)'
nhLabel   = r'$N_\mathrm{H}$ (cm$^{-2}$)'
densLabel = r'$n_\mathrm{H}$ (cm$^{-3}$)'
tempLabel = r'Temperature (K)'
fluxLabel = r'Local Flux (erg cm$^{-2}$ s$^{-1}$)'
abunLabel = r'$n(\mathrm{X})/n_\mathrm{H}$'
popLabel  = r'$n_i$ (cm$^{-3}$)'
coolLabel = r'$\Lambda$ (erg cm$^{-3}$ s$^{-1}$)'
heatLabel = r'$\Gamma$ (erg cm$^{-3}$ s$^{-1}$)'
rateLabel = r'Total Rate (erg cm$^{-3}$ s$^{-1}$)'
oprLabel  = r'$n($o-H$_2)/n($p-H$_2)$'

# Create the figure and specify the positions of the subplots
figure = pyplot.figure()
subplotPositions = gridspec.GridSpec(1, 1)

# Create the temperature plot
if temperaturePlot:
    print 'Creating temperature plot...'

    plotPosition = 1
    tickLabels = [True,True]
    axisLabels = [True,True]

    # Create the x- and y-axis data for each iteration
    xData = [] ; yData = [] ; dataLabel = []
    for iteration in range(len(temperatureArray)):
        if plotVersusAV:   xData.append([visualExtinction[i][rayID] for i in range(nParticles)])
        elif plotVersusNH: xData.append([totalColumn[i][rayID] for i in range(nParticles)])
        yData.append([temperatureArray[iteration][2],temperatureArray[iteration][3],temperatureArray[iteration][1]])

    dataLabel = [r'$T_\mathrm{min}$',r'$T_\mathrm{max}$',r'$T_\mathrm{new}$']
    color = ['#BBBBBB','#BBBBBB','black']

    # Specify the appropriate axis ranges and labels
    if plotVersusAV:
        xAxisRange = extinctionAxisRange
        xLabel = avLabel
    elif plotVersusNH:
        xAxisRange = columnAxisRange
        xLabel = nhLabel
    yAxisRange = temperatureAxisRange
    yLabel = tempLabel

    # Create the plot
    plottedLines, plot = create_plot(xData[0], yData[0], xAxisRange, yAxisRange, xLabel, yLabel, dataLabel, plotPosition, subplotPositions, logScale=[True,True], tickLabels=tickLabels, axisLabels=axisLabels, addLegend=False, color=color)
    line = plot.fill_between(xData[0], yData[0][0], yData[0][1], color=color[0])
    plottedLines.append(line)

# Create the density plot
if densityPlot:
    print 'Creating density plot...'

    if labelTickPositions.count(densityPlotPosition) > 0:
        tickLabels = [True,True]
    else:
        tickLabels = [False,False]
    if labelAxesPositions.count(densityPlotPosition) > 0:
        axisLabels = [True,True]
    else:
        axisLabels = [False,False]

    # Create the y-axis data
    xData = [] ; yData = [] ; dataLabel = []
    if plotVersusAV:
       xData = [visualExtinction[i][rayID] for i in range(nParticles)]
       xAxisRange = extinctionAxisRange
       xLabel = avLabel
    elif plotVersusNH:
       xData = [totalColumn[i][rayID] for i in range(nParticles)]
       xAxisRange = columnAxisRange
       xLabel = nhLabel
    yData.append(gasDensity)
    yData.append(dustDensity)
    yAxisRange = densityAxisRange
    yLabel = densLabel
    dataLabel = [r'$n_\mathrm{H}$',r'$n_\mathrm{grain}$']

    # Create the plot
    densityPlot = create_plot(xData, yData, xAxisRange, yAxisRange, xLabel, yLabel, dataLabel, densityPlotPosition, subplotPositions, logScale=[True,True], tickLabels=tickLabels, axisLabels=axisLabels, addLegend=False, color=['red','blue'])

# Create the FUV and X-ray flux plot
if fluxPlot:
    print 'Creating incident flux plot...'

    if labelTickPositions.count(fluxPlotPosition) > 0:
        tickLabels = [True,True]
    else:
        tickLabels = [False,False]
    if labelAxesPositions.count(fluxPlotPosition) > 0:
        axisLabels = [True,True]
    else:
        axisLabels = [False,False]

    # Create the y-axis data
    xData = [] ; yData = [] ; dataLabel = []
    if plotVersusAV:
       xData = [visualExtinction[i][rayID] for i in range(nParticles)]
       xAxisRange = extinctionAxisRange
       xLabel = avLabel
    elif plotVersusNH:
       xData = [totalColumn[i][rayID] for i in range(nParticles)]
       xAxisRange = columnAxisRange
       xLabel = nhLabel
    yData.append(fuvFlux)
    yData.append(xrayFlux)
    yAxisRange = fluxAxisRange
    yLabel = fluxLabel
    dataLabel = [r'$F_\mathrm{FUV}$',r'$F_\mathrm{X}$']

    # Create the plot
    fluxPlot = create_plot(xData, yData, xAxisRange, yAxisRange, xLabel, yLabel, dataLabel, fluxPlotPosition, subplotPositions, logScale=[True,True], tickLabels=tickLabels, axisLabels=axisLabels, addLegend=True, color=['red','blue'])

# Create the fractional abundance plot
if abundancePlot:
    print 'Creating abundance plot...'

    plotPosition = 1
    tickLabels = [True,True]
    axisLabels = [True,True]

    # Create the x- and y-axis data for each iteration
    xData = [] ; yData = [] ; dataLabel = []
    for iteration in range(len(abundanceArray)):
        if plotVersusAV:   xData.append([visualExtinction[i][rayID] for i in range(nParticles)])
        elif plotVersusNH: xData.append([totalColumn[i][rayID] for i in range(nParticles)])
        yData.append([])
        for species in speciesList:
            if speciesName.count(species) == 1:
                yData[iteration].append([abundanceArray[iteration][i][speciesName.index(species)] for i in range(nParticles)])

    # Create the data label for each species
    for species in speciesList:
        dataLabel.append(convert_species_name(species))

    # Specify the appropriate axis ranges and labels
    if plotVersusAV:
        xAxisRange = extinctionAxisRange
        xLabel = avLabel
    elif plotVersusNH:
        xAxisRange = columnAxisRange
        xLabel = nhLabel
    yAxisRange = abundanceAxisRange
    yLabel = abunLabel

    # Create the plot
    color = ['red','orange','green','blue','purple','black','red','orange','green','blue','purple','black']
    style = ['solid','solid','solid','solid','solid','solid','dashed','dashed','dashed','dashed','dashed','dashed']
    plottedLines, plot = create_plot(xData[0], yData[0], xAxisRange, yAxisRange, xLabel, yLabel, dataLabel, plotPosition, subplotPositions, logScale=[True,True], tickLabels=tickLabels, axisLabels=axisLabels, addLegend=True, color=color, style=style)

# Create the population density plot
if populationPlot:
    print 'Creating population density plot...'

    plotPosition = 1
    tickLabels = [True,True]
    axisLabels = [True,True]

    # Create the x- and y-axis data for each iteration
    # Replace zeroes with minimum values
    xData = [] ; yData = [] ; dataLabel = []
    for iteration in range(len(populationArray)):
        if plotVersusAV:   xData.append([visualExtinction[i][rayID] for i in range(nParticles)])
        elif plotVersusNH: xData.append([totalColumn[i][rayID] for i in range(nParticles)])
        yData.append([])
        for n in levelList:
            yData[iteration].append([populationArray[iteration][i][n] if populationArray[iteration][i][n] > 0 else populationAxisRange[0][0]*0.1 for i in range(nParticles)])

    # Create the data label for each level population
    for n in levelList:
        dataLabel.append(convert_species_name(populationName[n][2:populationName[n].index(',')])+'('+populationName[n][populationName[n].index(',')+1:])

    # Specify the appropriate axis ranges and labels
    if plotVersusAV:
        xAxisRange = extinctionAxisRange
        xLabel = avLabel
    elif plotVersusNH:
        xAxisRange = columnAxisRange
        xLabel = nhLabel
    yAxisRange = populationAxisRange
    yLabel = popLabel

    # Create the plot
    plottedLines, plot = create_plot(xData[0], yData[0], xAxisRange, yAxisRange, xLabel, yLabel, dataLabel, plotPosition, subplotPositions, logScale=[True,True], tickLabels=tickLabels, axisLabels=axisLabels, addLegend=True)

# Create the local emissivity plot
if emissivityPlot:
    print 'Creating emissivity plot...'

    plotPosition = 1
    tickLabels = [True,True]
    axisLabels = [True,True]

    # Create the x- and y-axis data for each iteration
    # Replace zeroes with minimum values
    xData = [] ; yData = [] ; dataLabel = []
    for iteration in range(len(emissivityArray)):
        if plotVersusAV:   xData.append([visualExtinction[i][rayID] for i in range(nParticles)])
        elif plotVersusNH: xData.append([totalColumn[i][rayID] for i in range(nParticles)])
        yData.append([])
        for n in lineList:
            yData[iteration].append([emissivityArray[iteration][i][n] if emissivityArray[iteration][i][n] > 0 else emissivityAxisRange[0][0]*0.1 for i in range(nParticles)])

    # Create the data label for each emission line
    for n in lineList:
       dataLabel.append(lineName[n])

    # Specify the appropriate axis ranges and labels
    if plotVersusAV:
        xAxisRange = extinctionAxisRange
        xLabel = avLabel
    elif plotVersusNH:
        xAxisRange = columnAxisRange
        xLabel = nhLabel
    yAxisRange = emissivityAxisRange
    yLabel = coolLabel

    # Create the plot
    plottedLines, plot = create_plot(xData[0], yData[0], xAxisRange, yAxisRange, xLabel, yLabel, dataLabel, plotPosition, subplotPositions, logScale=[True,True], tickLabels=tickLabels, axisLabels=axisLabels, addLegend=True)

# Create the cooling rate plot
if coolingPlot:
    print 'Creating cooling rate plot...'

    plotPosition = 1
    tickLabels = [True,True]
    axisLabels = [True,True]

    # Create the x- and y-axis data for each iteration
    # Replace zeroes with minimum values
    xData = [] ; yData = [] ; dataLabel = []
    for iteration in range(len(coolingArray)):
        if plotVersusAV:   xData.append([visualExtinction[i][rayID] for i in range(nParticles)])
        elif plotVersusNH: xData.append([totalColumn[i][rayID] for i in range(nParticles)])
        yData.append([])
        for n in range(nCool):
            yData[iteration].append([coolingArray[iteration][i][n] if coolingArray[iteration][i][n] > 0 else coolingAxisRange[0][0]*0.1 for i in range(nParticles)])

    # Specify the appropriate axis ranges and labels
    if plotVersusAV:
        xAxisRange = extinctionAxisRange
        xLabel = avLabel
    elif plotVersusNH:
        xAxisRange = columnAxisRange
        xLabel = nhLabel
    yAxisRange = coolingAxisRange
    yLabel = coolLabel
    dataLabel = coolantName

    # Create the plot
    plottedLines, plot = create_plot(xData[0], yData[0], xAxisRange, yAxisRange, xLabel, yLabel, dataLabel, plotPosition, subplotPositions, logScale=[True,True], tickLabels=tickLabels, axisLabels=axisLabels, addLegend=True)

# Create the heating rate plot
if heatingPlot:
    print 'Creating heating rates plot...'

    if labelTickPositions.count(heatingPlotPosition) > 0:
        tickLabels = [True,True]
    else:
        tickLabels = [False,False]
    if labelAxesPositions.count(heatingPlotPosition) > 0:
        axisLabels = [True,True]
    else:
        axisLabels = [False,False]

    # Create the y-axis data
    # Replace zeroes with minimum values
    xData = [] ; yData = [] ; dataLabel = []
    if plotVersusAV:
       xData = [visualExtinction[i][rayID] for i in range(nParticles)]
       xAxisRange = extinctionAxisRange
       xLabel = avLabel
    elif plotVersusNH:
       xData = [totalColumn[i][rayID] for i in range(nParticles)]
       xAxisRange = columnAxisRange
       xLabel = nhLabel
    for n in range(nHeat):
       yData.append([heatingRate[i][n] if heatingRate[i][n] > 0 else heatingAxisRange[0][0]*0.1 for i in range(nParticles)])
    yAxisRange = heatingAxisRange
    yLabel = heatLabel
    dataLabel = heatingName

    # Create the plot
    color = ['red','orange','green','blue','purple','red','orange','green','blue','purple','red','orange','green','blue','purple']
    style = ['solid','solid','solid','solid','solid','dashed','dashed','dashed','dashed','dashed','dotted','dotted','dotted','dotted','dotted']
    color[nHeat-1] = 'black'
    style[nHeat-1] = 'dotted'
    heatingPlot = create_plot(xData, yData, xAxisRange, yAxisRange, xLabel, yLabel, dataLabel, heatingPlotPosition, subplotPositions, logScale=[True,True], tickLabels=tickLabels, axisLabels=axisLabels, addLegend=False, color=color, style=style)
    if heatingPlotPosition == 2 or heatingPlotPosition == 4:
        pyplot.legend(bbox_to_anchor=(1.0,0.5), loc=6, frameon=False, prop={'size':10}, borderpad=0, labelspacing=0.4, handletextpad=0.1, borderaxespad=0.3, handlelength=2.4)

# Create the total heating/cooling rate plot
if heatCoolPlot:
    print 'Creating heating/cooling plot...'

    plotPosition = 1
    tickLabels = [True,True]
    axisLabels = [True,True]

    # Create the x- and y-axis data for each iteration
    xData = [] ; yData = [] ; dataLabel = []
    for iteration in range(len(heatCoolArray)):
        if plotVersusAV:   xData.append([visualExtinction[i][rayID] for i in range(nParticles)])
        elif plotVersusNH: xData.append([totalColumn[i][rayID] for i in range(nParticles)])
        yData.append(heatCoolArray[iteration])

    dataLabel = [r'$\Gamma$',r'$\Lambda$']
    color = ['red','blue']

    # Specify the appropriate axis ranges and labels
    if plotVersusAV:
        xAxisRange = extinctionAxisRange
        xLabel = avLabel
    elif plotVersusNH:
        xAxisRange = columnAxisRange
        xLabel = nhLabel
    yAxisRange = heatingAxisRange
    yLabel = rateLabel

    # Create the plot
    plottedLines, plot = create_plot(xData[0], yData[0], xAxisRange, yAxisRange, xLabel, yLabel, dataLabel, plotPosition, subplotPositions, logScale=[True,True], tickLabels=tickLabels, axisLabels=axisLabels, addLegend=True, color=color)

# Create the H2 ortho/para ratio plot
if opRatioPlot:
    print 'Creating H2 ortho/para plot...'

    if labelTickPositions.count(op_ratioPlotPosition) > 0:
        tickLabels = [True,True]
    else:
        tickLabels = [False,False]
    if labelAxesPositions.count(op_ratioPlotPosition) > 0:
        axisLabels = [True,True]
    else:
        axisLabels = [False,False]

    # Create the y-axis data
    xData = [] ; yData = [] ; dataLabel = []
    if plotVersusAV:
       xData = [visualExtinction[i][rayID] for i in range(nParticles)]
       xAxisRange = extinctionAxisRange
       xLabel = avLabel
    elif plotVersusNH:
       xData = [totalColumn[i][rayID] for i in range(nParticles)]
       xAxisRange = columnAxisRange
       xLabel = nhLabel
    yData.append(ortho_paraRatio)
    yAxisRange = op_ratioAxisRange
    yLabel = oprLabel
    dataLabel = [oprLabel]

    # Create the plot
    op_ratioPlot = create_plot(xData, yData, xAxisRange, yAxisRange, xLabel, yLabel, dataLabel, op_ratioPlotPosition, subplotPositions, logScale=[True,False], tickLabels=tickLabels, axisLabels=axisLabels, addLegend=False, color=['black'])

    # Create the y-axis data
    xData = [] ; yData = [] ; dataLabel = []
    if plotVersusAV:
       xData = [visualExtinction[i][rayID] for i in range(nParticles)]
       xAxisRange = extinctionAxisRange
       xLabel = avLabel
    elif plotVersusNH:
       xData = [totalColumn[i][rayID] for i in range(nParticles)]
       xAxisRange = columnAxisRange
       xLabel = nhLabel
    yData.append(paraDensity)
    yData.append(orthoDensity)
    yAxisRange = op_densityAxisRange
    yLabel = oprLabel
    dataLabel = [r'$n($p-H$_2)$',r'$n($o-H$_2)$']

#    # Create the plot
#    op_ratioPlot = create_plot(xData, yData, xAxisRange, yAxisRange, xLabel, yLabel, dataLabel, op_ratioPlotPosition, subplotPositions, logScale=[True,True], tickLabels=[False,False], axisLabels=[False,False], addLegend=True, color=['red','blue'])

# Add the iteration label to the plot and specify its appropriate format string
iterationFormat = 'Iteration %3d'
iterationLabel = pyplot.text(0.015, 0.02, '', color='black', size=12, horizontalalignment='left', verticalalignment='bottom', transform=plot.transAxes)

def init():
    if temperaturePlot:
        plottedLines.pop(-1).remove()
        line = plot.fill_between(xData[0], yData[0][0], yData[0][1], color=color[0])
        plottedLines.append(line)
    for n in range(len(yData[0])):
        plottedLines[n].set_data(xData[0], yData[0][n])
    iterationLabel.set_text('')
    return plottedLines, iterationLabel

def animate(i):
    if temperaturePlot:
        plottedLines.pop(-1).remove()
        line = plot.fill_between(xData[i], yData[i][0], yData[i][1], color=color[0])
        plottedLines.append(line)
    for n in range(len(yData[i])):
        plottedLines[n].set_data(xData[i], yData[i][n])
    iterationLabel.set_text(iterationFormat % (i))
    return plottedLines, iterationLabel

# Create the animation and save it as an MP4 file
print 'Creating animation...'
animation = animation.FuncAnimation(figure, animate, numpy.arange(0,len(yData)), interval=25, blit=False, init_func=init)
print 'Saving animation...'
animation.save('3D-PDR-Iterations.mp4', fps=10, codec='mpeg4', clear_temp=True)

stop = time.time()
duration = stop - start
print '\nFinished!'
print 'Run time =', duration, 'seconds'
