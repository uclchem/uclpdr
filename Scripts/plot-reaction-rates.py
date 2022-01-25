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
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FormatStrFormatter, LogLocator, MultipleLocator
from operator import mul as multiply

# Specify the path for the ASCII data files
path = ''

# Specify the prefix and suffix for the model output data files
prefix = 'rollig_V1'
suffix = '.out'

# Specify whether to plot visual extinction (A_V, True) or total column density (N_H, False) for the x-axis
plotVersusAV = True
plotVersusNH = not plotVersusAV

# Specify the species whose formation/destruction reactions are to be plotted
species = 'H+'

# Specify the depth index to use when sorting the rate coefficients
sortIndex = 0

# Specify the number of reaction rates to plot per figure
nLines = 5

# Specify the plot axis ranges and major/minor tickmark spacings
#extinctionAxisRange  = [[1E-7,1E+1], [100,100]]
extinctionAxisRange  = [[1E-4,2E+1], [10.,1.0]]
#columnAxisRange      = [[1E18,1E23], [10.,1.0]]
columnAxisRange      = [[1E20,1E25], [10.,1.0]]

# Specify which ray ID to plot against
rayID = 0

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

# Read the reaction rate coefficients
def read_rate_coefficients(inputfile):
    inputfile.seek(0)
    data = inputfile.readline().split()
    reactantNames = [reactionLabel.split(',') for reactionLabel in data[2:]]
    rateCoefficient = []
    data = inputfile.readline().split()
    while data != []:
        nReactions = len(data)-1
        rateCoefficient.append([number(data[i]) for i in range(1,len(data))])
        data = inputfile.readline().split()
    return nReactions, reactantNames, rateCoefficient

# Read the reactants and products for each reaction
def read_reactants_products(inputfile):
    inputfile.seek(0)
    reactantList = [] ; productList = []
    data = inputfile.readline().split(',')
    while len(data) != 0 and data != ['']:
        reactantList.append([species for species in data[1:4] if species != ''])
        productList.append([species for species in data[4:8] if species != ''])
        data = inputfile.readline().split(',')
    return reactantList, productList

# Calculate the reaction rates (in s^-1) from the rate coefficients and fractional abundances
def calculate_reaction_rates(reactantList, speciesName, rateCoefficient, fractionalAbundance, gasDensity):
    reactionRate = rateCoefficient
    for j in range(len(reactantList)):
        abundanceIndex = [speciesName.index(reactant) for reactant in reactantList[j] if speciesName.count(reactant) == 1]
        for i in range(len(rateCoefficient)):
            reactionRate[i][j] = rateCoefficient[i][j]*(gasDensity[i]**(len(abundanceIndex)-1))*reduce(multiply, [fractionalAbundance[i][abundanceIndex[k]] for k in range(len(abundanceIndex))])
    return reactionRate

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

    # Deal with special reactant names
    if name == '#':      return r'\#'
    if name == 'e-':     return r'$e^{-}$'
    if name == 'PHOTON': return r'$h\nu$'
    if name == 'CRPHOT': return r'$\gamma_{_\mathrm{CR}}$'
    if name == 'CRP':    return r'\textsc{crp}'
    if name == 'XRAY':   return r'$\gamma_{_\mathrm{XR}}$'
    if name == 'XRSEC':  return r'$e^{-}_{_\mathrm{XR}}$'
    if name == 'XRLYA':  return r'$\gamma_{_\mathrm{XR,\,Ly\alpha}}$'
    if name == 'XRPHOT': return r'$\gamma_{_\mathrm{XR,\,LyW}}$'

    string = ''
    for i in range(len(name)):

        # Handle the special case of ortho/para states
        if name[i:i+2] == 'p-': string += r'p-' ; continue
        if name[i:i+2] == 'o-': string += r'o-' ; continue
        if name[i-1:i+1] == 'p-': continue
        if name[i-1:i+1] == 'o-': continue

        # Handle the special case of doubly ionised species
        if name[i:i+2] == '++': string += r'$^{++}$' ; continue
        if name[i:i+2] == '--': string += r'$^{--}$' ; continue
        if name[i-1:i+1] == '++': continue
        if name[i-1:i+1] == '--': continue

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
        plot = pyplot.plot(xData, yData[n], color=color[n], linestyle=style[n], linewidth=1.0, label=dataLabel[n])

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
        pyplot.legend(mode='expand', loc='best', frameon=False, prop={'size':10}, borderpad=0, labelspacing=0.1, handletextpad=0.1, borderaxespad=0.6, handlelength=2.4)

    return plot

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

print 'Reading rate coefficients...'
dataFilename = prefix+'.reac.out'
input = codecs.open(path+dataFilename, encoding='utf-8', mode='r')
nReactions, reactantNames, rateCoefficient = read_rate_coefficients(input)
input.close()

print 'Reading fractional abundances...'
dataFilename = prefix+'.abun'+suffix
input = codecs.open(path+dataFilename, encoding='utf-8', mode='r')
nSpecies, speciesName, fractionalAbundance = read_fractional_abundances(input)
input.close()

print 'Reading reactants and products...'
dataFilename = prefix+'.rates.inp'
input = codecs.open(path+dataFilename, encoding='utf-8', mode='r')
reactantList, productList = read_reactants_products(input)
input.close()

# Open the output PDF file to store the figures
pdfFile = PdfPages('Reaction-Rates.pdf')

# Convert the visual extinction values (A_V) to total column densities (N_H)
totalColumn = [[A_V_ray/6.289E-22 for A_V_ray in A_V] for A_V in visualExtinction]

# Calculate the reaction rates
print 'Calculating reaction rates...'
reactionRate = calculate_reaction_rates(reactantList, speciesName, rateCoefficient, fractionalAbundance, gasDensity)

# Create sorted reaction lists for the formation and destruction of a desired species
print 'Finding formation/destruction reactions...'
formationIndices = [] ; destructionIndices = [] ; reactionIndices = []
for n in range(nReactions):
    if species in productList[n]:
        formationIndices.append((reactionRate[sortIndex][n], n))
        productList[n].pop(productList[n].index(species))
        productList[n].insert(0, species)
    if species in reactantList[n]:
        destructionIndices.append((reactionRate[sortIndex][n], n))
        reactantList[n].pop(reactantList[n].index(species))
        reactantList[n].insert(0, species)
formationIndices.sort(reverse=True) ; destructionIndices.sort(reverse=True)
reactionIndices.extend([formationIndices[i][1] for i in range(nLines*(len(formationIndices)//nLines))])
reactionIndices.extend([destructionIndices[i][1] for i in range(len(destructionIndices))])
nReactions = len(reactionIndices)

# Specify the axis labels
avLabel   = r'$A_V$ (mag)'
nhLabel   = r'$N_\mathrm{H}$ (cm$^{-2}$)'
rateLabel = r'$R_i$ (s$^{-1}$)'

# Determine the total number of figures needed to plot all reaction rates
nPlots = int(math.ceil(float(nReactions)/float(nLines)))

# Loop over all reactions, plotting nLines per figure and four figures per page
for plotNumber in range(nPlots):

    # Clear and set up a new page for each four figures
    if plotNumber % 4 == 0:

        # Add the previous figures to the output PDF file
        if plotNumber > 0:
            print 'Saving current page to PDF file...'
            pdfFile.savefig()

        # Close any previously open plot(s)
        pyplot.close()

        # Create the figure and specify the positions of the subplots
        figure = pyplot.figure(1)
        subplotPositions = gridspec.GridSpec(2, 2)
        subplotPositions.update(left=0.15, right=0.85, bottom=0.08, top=0.98, wspace=0.33, hspace=0.25)
        labelTickPositions = range(1,5)
        labelAxesPositions = range(1,5)

    # Create the reaction rate plot
    print 'Creating reaction rate plot (' + str(plotNumber+1) + '/' + str(nPlots) + ')...'

    plotPosition = (plotNumber % 4) + 1

    if labelTickPositions.count(plotPosition) > 0:
        tickLabels = [True,True]
    else:
        tickLabels = [False,False]

    if labelAxesPositions.count(plotPosition) > 0:
        axisLabels = [True,True]
    else:
        axisLabels = [False,False]

    # Create the x- and y-axis data
    xData = [] ; yData = [] ; dataLabel = []
    if plotVersusAV:
        xData = [visualExtinction[i][rayID] for i in range(nParticles)]
        xAxisRange = extinctionAxisRange
        xLabel = avLabel
    elif plotVersusNH:
        xData = [totalColumn[i][rayID] for i in range(nParticles)]
        xAxisRange = columnAxisRange
        xLabel = nhLabel
    for reactionNumber in range(plotNumber*nLines,min((plotNumber+1)*nLines, nReactions)):
        index = reactionIndices[reactionNumber]
        yData.append([max(reactionRate[i][index], 1.0E-99) for i in range(nParticles)])
        dataLabel.append(r'\makebox[20mm][l]{' + r' $+$ '.join([convert_species_name(reactant) for reactant in reactantList[index]]) + r'}' \
                       + r' $\rightarrow$ ' \
                       + r'\makebox[20mm][l]{' + r' $+$ '.join([convert_species_name(product) for product in productList[index]]) + r'}')
    yLabel = rateLabel

    # Determine the appropriate y-axis range
#     yAxisMax = math.ceil(math.log10(max([value for dataset in yData for value in dataset])))
    yAxisMax = math.ceil(math.log10(max([dataset[sortIndex] for dataset in yData])))
    yAxisMin = math.floor(math.log10(min([value for dataset in yData for value in dataset])))
    if yAxisMax < -30: yAxisMax = -30
    if yAxisMax - yAxisMin > 0:
        yAxisRange = [[10**(yAxisMin),10**(yAxisMax+2)], [10,1]]
    if yAxisMax - yAxisMin < 3:
        yAxisRange = [[10**(yAxisMax-2),10**(yAxisMax+2)], [10,1]]
    if yAxisMax - yAxisMin > 7:
        yAxisRange = [[10**(yAxisMax-8),10**(yAxisMax+2)], [100,100]]

    # Create the plot
    if len(dataLabel) % 6 == 0:
        color = ['red','orange','green','blue','purple','black','red','orange','green','blue','purple','black']
        style = ['solid','solid','solid','solid','solid','dashed','dashed','dashed','dashed','dashed','dashed']
    else:
        color = ['red','orange','green','blue','purple','red','orange','green','blue','purple','red','orange','green','blue','purple']
        style = ['solid','solid','solid','solid','solid','dashed','dashed','dashed','dashed','dashed','dotted','dotted','dotted','dotted','dotted']
    reacRatePlot = create_plot(xData, yData, xAxisRange, yAxisRange, xLabel, yLabel, dataLabel, plotPosition, subplotPositions, logScale=[True,True], tickLabels=tickLabels, axisLabels=axisLabels, addLegend=True, color=color, style=style)

    # Set the background to be transparent if the figure is to be used in a presentation
    if useForPresentation:
        figure.patch.set_alpha(0)

# Add the final figures to the output PDF file
print 'Saving final page to PDF file...'
pdfFile.savefig()

# Close the output PDF file
pdfFile.close()
pyplot.close()

stop = time.time()
duration = stop - start
print '\nFinished!'
print 'Run time =', duration, 'seconds'
