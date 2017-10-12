#! /usr/bin/python
#
import codecs
import math
import os
import time
import numpy
import matplotlib
matplotlib.use('pdf')
import matplotlib.colors as colors
import matplotlib.pyplot as pyplot
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter, LogLocator, MultipleLocator

# Specify the path for the ASCII data files
path = ''

# Specify the prefix and suffix for the model output data files
prefix = 'rollig_V1'
suffix = '.out'

# Specify the location for each plot (-1 to exclude a plot from the figure)
fluxPlotPosition        = -1
densityPlotPosition     = -1
densTempPlotPosition    = -1
temperaturePlotPosition = 1
abundancePlotPosition   = 2
populationPlotPosition  = -1
emissivityPlotPosition  = -1
opacityPlotPosition     = -1
coolingPlotPosition     = 3
heatingPlotPosition     = 4
heatCoolPlotPosition    = -1
intensityPlotPosition   = -1
opRatioPlotPosition     = -1

# Specify whether to plot visual extinction (A_V, True) or total column density (N_H, False) for the x-axis
plotVersusAV = True
plotVersusNH = not plotVersusAV

# Specify if the temperatures are to be plotted on a logarithmic scale
plotLogTemperature = True

# Specify the plot axis ranges and major/minor tickmark spacings
extinctionAxisRange  = [[1E-4,2E+1], [10.,1.0]]
columnAxisRange      = [[1E20,1E25], [10.,1.0]]
fluxAxisRange        = [[1E-5,1E+3], [10.,1.0]]
densityAxisRange     = [[1E+1,1E+6], [10.,1.0]]
temperatureAxisRange = [[1E+0,1E+3], [10.,1.0]]
abundanceAxisRange   = [[1E-7,2E+1], [10.,1.0]]
populationAxisRange  = [[1E-10,1E+0], [100,100]]
emissivityAxisRange  = [[1E-28,1E-18], [100,100]]
coolingAxisRange     = [[1E-28,1E-18], [100,100]]
heatingAxisRange     = [[1E-28,1E-18], [100,100]]
opacityAxisRange     = [[1E-3,1E+2], [10.,1.0]]
intensityAxisRange   = [[1E-8,1E-3], [10.,1.0]]
opDensityAxisRange   = [[1E-5,1E+5], [100,100]]
opRatioAxisRange     = [[+0.0,+6.0], [100,100]]

if not plotLogTemperature:
    temperatureAxisRange = [[0.0,100], [20.,5.0]]

# Specify which ray ID to plot against
rayID = 0

# Specify the species whose abundances should be plotted
speciesList = ['H','H2','C+','C','CO','e-'] # PDR
# speciesList = ['H','H2','O','C','CO','e-','H+','O+','C+'] # XDR

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
        lineName[i] = lineName[i].replace('III','{\sc\,iii}')
        lineName[i] = lineName[i].replace('II','{\sc\,ii}')
        lineName[i] = lineName[i].replace('I','{\sc\,i}')
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

# Read the line opacities
def read_line_opacities(inputfile):
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
        lineName[i] = lineName[i].replace('III','{\sc\,iii}')
        lineName[i] = lineName[i].replace('II','{\sc\,ii}')
        lineName[i] = lineName[i].replace('I','{\sc\,i}')
        lineName[i] = lineName[i].replace('um','$\upmu$m')
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
        coolantName[i] = coolantName[i].replace('III','{\sc\,iii}')
        coolantName[i] = coolantName[i].replace('II','{\sc\,ii}')
        coolantName[i] = coolantName[i].replace('I','{\sc\,i}')
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
        if style[n] in ['solid','dashed','dotted','dashdot']:
            plot, = pyplot.plot(xData, yData[n], color=color[n], linestyle=style[n], linewidth=1.0, label=dataLabel[n])
        else:
            plot, = pyplot.plot(xData, yData[n], color='#DDDDDD', linestyle='solid', linewidth=1.0, label=dataLabel[n],
                                marker=style[n], markeredgewidth=0.6, markersize=10, markerfacecolor=color[n])
        plottedLines.append(plot)

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

    # Add a legend, if requested
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

    # Set the tickmarks, axes and labels to white if the figure is to be used in a presentation
    if useForPresentation:
        pyplot.setp(subplot.xaxis.get_majorticklines(), color='white')
        pyplot.setp(subplot.xaxis.get_minorticklines(), color='white')
        pyplot.setp(subplot.yaxis.get_majorticklines(), color='white')
        pyplot.setp(subplot.yaxis.get_minorticklines(), color='white')
        pyplot.setp(subplot.spines.values(), color='white')
        subplot.xaxis.set_tick_params(labelcolor='white')
        subplot.yaxis.set_tick_params(labelcolor='white')
        if axisLabels[0]:
            subplot.xaxis.label.set_color('white')
        if axisLabels[1]:
            subplot.yaxis.label.set_color('white')
        if addLegend:
            pyplot.setp(subplot.get_legend().get_texts(), color='white')
    return subplot, plottedLines

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

print 'Reading fractional abundances...'
dataFilename = prefix+'.abun'+suffix
input = codecs.open(path+dataFilename, encoding='utf-8', mode='r')
nSpecies, speciesName, fractionalAbundance = read_fractional_abundances(input)
input.close()

print 'Reading population densities...'
dataFilename = prefix+'.pop'+suffix
input = codecs.open(path+dataFilename, encoding='utf-8', mode='r')
nLevels, populationName, populationDensity = read_population_densities(input)
input.close()

print 'Reading line emissivities...'
dataFilename = prefix+'.emis'+suffix
input = codecs.open(path+dataFilename, encoding='utf-8', mode='r')
nLines, lineName, localEmissivity = read_local_emissivities(input)
input.close()

print 'Reading line opacities...'
dataFilename = prefix+'.tau_'+'0'*(int(math.log10(nRays))-int(math.log10(rayID+1)))+str(rayID+1)+suffix
input = codecs.open(path+dataFilename, encoding='utf-8', mode='r')
nLines, lineName, lineOpacity = read_line_opacities(input)
input.close()

print 'Reading cooling rates...'
dataFilename = prefix+'.cool'+suffix
input = codecs.open(path+dataFilename, encoding='utf-8', mode='r')
nCool, coolantName, coolingRate = read_cooling_rates(input)
input.close()

print 'Reading heating rates...'
dataFilename = prefix+'.heat'+suffix
input = codecs.open(path+dataFilename, encoding='utf-8', mode='r')
nHeat, heatingName, heatingRate = read_heating_rates(input)
input.close()

# Convert the visual extinction values (A_V) to total column densities (N_H)
totalColumn = [[A_V_ray/6.289E-22 for A_V_ray in A_V] for A_V in visualExtinction]

# Convert the FUV values from multiples of the Draine field to fluxes in erg cm^-2 s^-1
fuvFlux = [chi*2.74E-3 for chi in fuvFlux]

# Add gas-grain collisions to the list of coolants when it acts as a negative heating rate
index = heatingName.index('Gas-Grain')
if min(heatingRate[:][index]) < 0:
    nCool += 1
    coolantName.insert(-1,r'Grain')
    for n in range(nParticles):
        coolingRate[n].insert(-1,-heatingRate[n][index])
        coolingRate[n][-1] = coolingRate[n][-1] + max(coolingRate[n][-2], 0)
        heatingRate[n][-1] = heatingRate[n][-1] + max(coolingRate[n][-2], 0)

# Calculate the integrated intensities of the main PDR diagnostic lines
# Assume the emerging lines emit over a full hemisphere (2pi steradians)
emissionLines = [] ; integratedIntensity = []
for lineLabel in emissionList:
    for n in range(nLines):
        if lineLabel in lineName[n]:
            emissionLines.append(lineName[n].replace(' ',' \n ').replace('(',' \n ('))
            integratedIntensity.append(sum([0.5*(localEmissivity[i+1][n]+localEmissivity[i][n])*(particleCoordinates[i+1][0]-particleCoordinates[i][0]) for i in range(nParticles-1)])/(2*math.pi))

# Calculate the ortho-to-para ratio of H2
paraDensity = [] ; orthoDensity = []
for n in range(nParticles):
    paraDensity.append(sum([populationDensity[n][i] for i in range(nLevels) if populationName[i].count('p-H2')]))
    orthoDensity.append(sum([populationDensity[n][i] for i in range(nLevels) if populationName[i].count('o-H2')]))
if min(paraDensity) > 0:
    orthoParaRatio = [orthoDensity[i]/paraDensity[i] for i in range(nParticles)]
else:
    orthoParaRatio = [3.0 for i in range(nParticles)]

# Specify the axis labels
avLabel   = r'$A_V$ (mag)'
nhLabel   = r'$N_\mathrm{H}$ (cm$^{-2}$)'
fluxLabel = r'Flux (erg cm$^{-2}$ s$^{-1}$)'
densLabel = r'$n_\mathrm{H}$ (cm$^{-3}$)'
tempLabel = r'Temperature (K)'
abunLabel = r'$n(\mathrm{X})/n_\mathrm{H}$'
coolLabel = r'$\Lambda$ (erg cm$^{-3}$ s$^{-1}$)'
heatLabel = r'$\Gamma$ (erg cm$^{-3}$ s$^{-1}$)'
lineLabel = r'$I$ (erg s$^{-1}$ cm$^{-2}$ sr$^{-1}$)'
popLabel  = r'$n_i$ (cm$^{-3}$)'
tauLabel  = r'$\tau$'
oprLabel  = r'$n($o-H$_2)/n($p-H$_2)$'

# Create the figure and specify the positions of the subplots
figure = pyplot.figure(1)
subplotPositions = gridspec.GridSpec(2, 2)
subplotPositions.update(left=0.15, right=0.85, bottom=0.08, top=0.98, wspace=0.33, hspace=0.25)
labelTickPositions = range(1,5)
labelAxesPositions = range(1,5)

print

# Create the FUV and X-ray flux plot
if fluxPlotPosition > 0:
    print 'Creating local flux plot...'

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
    subplot, fluxPlot = create_plot(xData, yData, xAxisRange, yAxisRange, xLabel, yLabel, dataLabel, fluxPlotPosition, subplotPositions, logScale=[True,True], tickLabels=tickLabels, axisLabels=axisLabels, addLegend=True, color=['red','blue'])

    # Set the background to be transparent if the figure is to be used in a presentation
    if useForPresentation:
        subplot.patch.set_alpha(0)

# Create the density plot
if densityPlotPosition > 0:
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
    subplot, densityPlot = create_plot(xData, yData, xAxisRange, yAxisRange, xLabel, yLabel, dataLabel, densityPlotPosition, subplotPositions, logScale=[True,True], tickLabels=tickLabels, axisLabels=axisLabels, addLegend=False, color=['red','blue'])

    # Set the background to be transparent if the figure is to be used in a presentation
    if useForPresentation:
        subplot.patch.set_alpha(0)

# Create the temperature plot
if temperaturePlotPosition > 0:
    print 'Creating temperature plot...'

    if labelTickPositions.count(temperaturePlotPosition) > 0:
        tickLabels = [True,True]
    else:
        tickLabels = [False,False]
    if labelAxesPositions.count(temperaturePlotPosition) > 0:
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
    yData.append(gasTemperature)
    yData.append(dustTemperature)
    yAxisRange = temperatureAxisRange
    yLabel = tempLabel
    dataLabel = [r'$T_\mathrm{gas}$',r'$T_\mathrm{dust}$']

    # Create the plot
    subplot, temperaturePlot = create_plot(xData, yData, xAxisRange, yAxisRange, xLabel, yLabel, dataLabel, temperaturePlotPosition, subplotPositions, logScale=[True,plotLogTemperature], tickLabels=tickLabels, axisLabels=axisLabels, addLegend=True, color=['red','blue'])

    # Set the background to be transparent if the figure is to be used in a presentation
    if useForPresentation:
        subplot.patch.set_alpha(0)

# Create the density-temperature combo plot
if densTempPlotPosition > 0:
    print 'Creating density and temperature plot...'

    if labelTickPositions.count(densTempPlotPosition) > 0:
        tickLabels = [True,True]
    else:
        tickLabels = [False,False]
    if labelAxesPositions.count(densTempPlotPosition) > 0:
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
    yData.append(gasTemperature)
    yData.append(dustTemperature)
    yData.append([0]*nParticles)
    yAxisRange = temperatureAxisRange
    yLabel = tempLabel
    dataLabel = [r'$T_\mathrm{gas}$',r'$T_\mathrm{dust}$',r'$n_\mathrm{H}$']

    # Create the plot
    subplot, densTempPlot = create_plot(xData, yData, xAxisRange, yAxisRange, xLabel, yLabel, dataLabel, densTempPlotPosition, subplotPositions, logScale=[True,plotLogTemperature], tickLabels=tickLabels, axisLabels=axisLabels, addLegend=True, color=['red','blue','black'], style=['solid','solid','dashed'])

    # Add the density profile with a separate y-axis scale on the right-hand side of the plot
    yData[-1] = gasDensity
    yAxisRange = densityAxisRange
    yLabel = r'Density (cm$^{-3}$)'
    rhsAxis = subplot.twinx()
    rhsAxis.plot(xData, yData[2], color='black', linestyle='dashed', linewidth=1.0, label=dataLabel[2])

    # Set the axes labels, scales, ranges and tickmarks
    rhsAxis.set_ylabel(yLabel)
    rhsAxis.set_xscale('log')
    rhsAxis.axis([xAxisRange[0][0], xAxisRange[0][1], yAxisRange[0][0], yAxisRange[0][1]])
    rhsAxis.xaxis.set_major_locator(LogLocator(base=xAxisRange[1][0]))
    rhsAxis.xaxis.set_minor_locator(LogLocator(base=xAxisRange[1][0], subs=range(int(xAxisRange[1][0]/xAxisRange[1][1]))))
    rhsAxis.set_yscale('log')
    rhsAxis.yaxis.set_major_locator(LogLocator(base=yAxisRange[1][0]))
    rhsAxis.yaxis.set_minor_locator(LogLocator(base=yAxisRange[1][0], subs=range(int(yAxisRange[1][0]/yAxisRange[1][1]))))

#    # Plot the density profile on a linear y-axis scale
#     rhsAxis.yaxis.set_major_locator(MultipleLocator(yAxisRange[1][0]))
#     rhsAxis.yaxis.set_minor_locator(MultipleLocator(yAxisRange[1][1]))
#     tickValues = rhsAxis.yaxis.get_majorticklocs()
#     tickLabels = [r'%i$\!\times\!$$10^{%i}$'% (value/(10**int(math.log10(tickValues[1]))), int(math.log10(tickValues[1]))) for value in tickValues]
#     rhsAxis.yaxis.set_ticklabels(tickLabels)

    # Set the background to be transparent if the figure is to be used in a presentation
    if useForPresentation:
        pyplot.setp(rhsAxis.yaxis.get_majorticklines(), color='white')
        pyplot.setp(rhsAxis.yaxis.get_minorticklines(), color='white')
        rhsAxis.yaxis.set_tick_params(labelcolor='white')
        rhsAxis.yaxis.label.set_color('white')
        subplot.patch.set_alpha(0)

# Create the fractional abundance plot
if abundancePlotPosition > 0:
    print 'Creating abundance plot...'

    if labelTickPositions.count(abundancePlotPosition) > 0:
        tickLabels = [True,True]
    else:
        tickLabels = [False,False]
    if labelAxesPositions.count(abundancePlotPosition) > 0:
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
    for species in speciesList:
        if speciesName.count(species) == 1:
            yData.append([fractionalAbundance[i][speciesName.index(species)] for i in range(nParticles)])
            dataLabel.append(convert_species_name(species))
    yAxisRange = abundanceAxisRange
    yLabel = abunLabel

    # Create the plot
    if len(dataLabel)%6 == 0:
        color = ['red','orange','green','blue','purple','black','red','orange','green','blue','purple','black']
        style = ['solid','solid','solid','solid','solid','dashed','dashed','dashed','dashed','dashed','dashed']
    else:
        color = ['red','orange','green','blue','purple','red','orange','green','blue','purple','red','orange','green','blue','purple']
        style = ['solid','solid','solid','solid','solid','dashed','dashed','dashed','dashed','dashed','dotted','dotted','dotted','dotted','dotted']
    subplot, abundancePlot = create_plot(xData, yData, xAxisRange, yAxisRange, xLabel, yLabel, dataLabel, abundancePlotPosition, subplotPositions, logScale=[True,True], tickLabels=tickLabels, axisLabels=axisLabels, addLegend=True, color=color, style=style)

    # Set the background to be transparent if the figure is to be used in a presentation
    if useForPresentation:
        subplot.patch.set_alpha(0)

# Create the local emissivity plot
if emissivityPlotPosition > 0:
    print 'Creating emissivity plot...'

    if labelTickPositions.count(emissivityPlotPosition) > 0:
        tickLabels = [True,True]
    else:
        tickLabels = [False,False]
    if labelAxesPositions.count(emissivityPlotPosition) > 0:
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
    for n in lineList:
        yData.append([localEmissivity[i][n] if localEmissivity[i][n] > 0 else emissivityAxisRange[0][0]*0.1 for i in range(nParticles)])
        dataLabel.append(lineName[n])
    yAxisRange = emissivityAxisRange
    yLabel = coolLabel

    # Create the plot
    subplot, emissivityPlot = create_plot(xData, yData, xAxisRange, yAxisRange, xLabel, yLabel, dataLabel, emissivityPlotPosition, subplotPositions, logScale=[True,True], tickLabels=tickLabels, axisLabels=axisLabels, addLegend=True)

    # Set the background to be transparent if the figure is to be used in a presentation
    if useForPresentation:
        subplot.patch.set_alpha(0)

# Create the line opacity plot
if opacityPlotPosition > 0:
    print 'Creating line opacity plot...'

    if labelTickPositions.count(opacityPlotPosition) > 0:
        tickLabels = [True,True]
    else:
        tickLabels = [False,False]
    if labelAxesPositions.count(opacityPlotPosition) > 0:
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
    for n in lineList:
        yData.append([lineOpacity[i][n] if lineOpacity[i][n] > 0 else opacityAxisRange[0][0]*0.1 for i in range(nParticles)])
        dataLabel.append(lineName[n])
    yAxisRange = opacityAxisRange
    yLabel = tauLabel

    # Create the plot
    subplot, opacityPlot = create_plot(xData, yData, xAxisRange, yAxisRange, xLabel, yLabel, dataLabel, opacityPlotPosition, subplotPositions, logScale=[True,True], tickLabels=tickLabels, axisLabels=axisLabels, addLegend=True)

    # Set the background to be transparent if the figure is to be used in a presentation
    if useForPresentation:
        subplot.patch.set_alpha(0)

# Create the cooling rate plot
if coolingPlotPosition > 0:
    print 'Creating cooling rate plot...'

    if labelTickPositions.count(coolingPlotPosition) > 0:
        tickLabels = [True,True]
    else:
        tickLabels = [False,False]
    if labelAxesPositions.count(coolingPlotPosition) > 0:
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
    for n in range(nCool):
        yData.append([coolingRate[i][n] if coolingRate[i][n] > 0 else coolingAxisRange[0][0]*0.1 for i in range(nParticles)])
    yAxisRange = coolingAxisRange
    yLabel = coolLabel
    dataLabel = coolantName

    # Create the plot
    color = ['red','orange','green','blue','purple','red','orange','green','blue','purple','red','orange','green','blue','purple']
    style = ['solid','solid','solid','solid','solid','dashed','dashed','dashed','dashed','dashed','dotted','dotted','dotted','dotted','dotted']
    color[nCool-1] = 'black'
    style[nCool-1] = 'dashed'
    subplot, coolingPlot = create_plot(xData, yData, xAxisRange, yAxisRange, xLabel, yLabel, dataLabel, coolingPlotPosition, subplotPositions, logScale=[True,True], tickLabels=tickLabels, axisLabels=axisLabels, addLegend=True, color=color, style=style)

    # Set the background to be transparent if the figure is to be used in a presentation
    if useForPresentation:
        subplot.patch.set_alpha(0)

# Create the heating rate plot
if heatingPlotPosition > 0:
    print 'Creating heating rate plot...'

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
    style[nHeat-1] = 'dashed'
    subplot, heatingPlot = create_plot(xData, yData, xAxisRange, yAxisRange, xLabel, yLabel, dataLabel, heatingPlotPosition, subplotPositions, logScale=[True,True], tickLabels=tickLabels, axisLabels=axisLabels, addLegend=False, color=color, style=style)

    # Add a legend on the right-hand side of the plot
    if heatingPlotPosition == 2 or heatingPlotPosition == 4:
        pyplot.legend(bbox_to_anchor=(1.0,0.5), loc=6, frameon=False, prop={'size':10}, borderpad=0, labelspacing=0.4, handletextpad=0.1, borderaxespad=0.3, handlelength=2.4)

    # Set the background to be transparent if the figure is to be used in a presentation
    if useForPresentation:
        subplot.patch.set_alpha(0)

# Create the heating/cooling comparison plot
if heatCoolPlotPosition > 0:
    print 'Creating heating/cooling plot...'

    if labelTickPositions.count(heatCoolPlotPosition) > 0:
        tickLabels = [True,True]
    else:
        tickLabels = [False,False]
    if labelAxesPositions.count(heatCoolPlotPosition) > 0:
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
    yData.append([heatingRate[i][-1] if heatingRate[i][-1] > 0 else heatingAxisRange[0][0]*0.1 for i in range(nParticles)])
    yData.append([coolingRate[i][-1] if coolingRate[i][-1] > 0 else heatingAxisRange[0][0]*0.1 for i in range(nParticles)])
    yAxisRange = heatingAxisRange
    yLabel = heatLabel
    dataLabel = ['Heating','Cooling']

    # Create the plot
    color = ['red','blue']
    style = ['solid','solid']
    subplot, heatCoolPlot = create_plot(xData, yData, xAxisRange, yAxisRange, xLabel, yLabel, dataLabel, heatCoolPlotPosition, subplotPositions, logScale=[True,True], tickLabels=tickLabels, axisLabels=axisLabels, addLegend=True, color=color, style=style)

    # Set the background to be transparent if the figure is to be used in a presentation
    if useForPresentation:
        subplot.patch.set_alpha(0)

# Create the integrated intensity plot
if intensityPlotPosition > 0:
    print 'Creating intensity plot...'

    if labelTickPositions.count(intensityPlotPosition) > 0:
        tickLabels = [True,True]
    else:
        tickLabels = [False,False]
    if labelAxesPositions.count(intensityPlotPosition) > 0:
        axisLabels = [True,True]
    else:
        axisLabels = [False,False]

    # Create the y-axis data
    # Replace zeroes with minimum values
    xData = range(1,len(emissionLines)+1)
    yData = [integratedIntensity]
    xAxisRange = [[0.5,len(emissionLines)+0.5], [1,1]]
    yAxisRange = intensityAxisRange
    xLabel = ''
    yLabel = lineLabel
    dataLabel = ['']*len(emissionLines)

    # Create the plot
    subplot, intensityPlot, = create_plot(xData, yData, xAxisRange, yAxisRange, xLabel, yLabel, dataLabel, intensityPlotPosition, subplotPositions, logScale=[False,True], tickLabels=tickLabels, axisLabels=axisLabels, addLegend=False, color=['#44AAFF'], style=['*'])
    subplot.xaxis.set_ticklabels(emissionLines, size=9)

    # Set the background to be transparent if the figure is to be used in a presentation
    if useForPresentation:
        subplot.patch.set_alpha(0)

    # Print out the values for common diagnostic line ratios
    print ' log(63um/146um) =%+7.3f' % math.log10(integratedIntensity[0]/integratedIntensity[1])
    print ' log(158um/63um) =%+7.3f' % math.log10(integratedIntensity[2]/integratedIntensity[0])

# Create the H2 ortho/para ratio plot
if opRatioPlotPosition > 0:
    print 'Creating H2 ortho/para plot...'

    if labelTickPositions.count(opRatioPlotPosition) > 0:
        tickLabels = [True,True]
    else:
        tickLabels = [False,False]
    if labelAxesPositions.count(opRatioPlotPosition) > 0:
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
    yData.append(orthoParaRatio)
    yAxisRange = opRatioAxisRange
    yLabel = oprLabel
    dataLabel = [oprLabel]

    # Create the plot
    subplot, opRatioPlot = create_plot(xData, yData, xAxisRange, yAxisRange, xLabel, yLabel, dataLabel, opRatioPlotPosition, subplotPositions, logScale=[True,False], tickLabels=tickLabels, axisLabels=axisLabels, addLegend=False, color=['black'])

    # Set the background to be transparent if the figure is to be used in a presentation
    if useForPresentation:
        subplot.patch.set_alpha(0)

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
    yAxisRange = opDensityAxisRange
    yLabel = oprLabel
    dataLabel = [r'$n($p-H$_2)$',r'$n($o-H$_2)$']

#     # Create the plot
#     subplot, opRatioPlot = create_plot(xData, yData, xAxisRange, yAxisRange, xLabel, yLabel, dataLabel, opRatioPlotPosition, subplotPositions, logScale=[True,True], tickLabels=[False,False], axisLabels=[False,False], addLegend=True, color=['red','blue'])
#
#     # Set the background to be transparent if the figure is to be used in a presentation
#     if useForPresentation:
#         subplot.patch.set_alpha(0)

# Set the background to be transparent if the figure is to be used in a presentation
if useForPresentation:
    figure.patch.set_alpha(0)

# Save the figure to the output PDF file
print 'Saving figure to PDF file...'
figure.savefig('Model-Results.pdf', transparent=useForPresentation)
pyplot.close()

stop = time.time()
duration = stop - start
print '\nFinished!'
print 'Run time =', duration, 'seconds'
