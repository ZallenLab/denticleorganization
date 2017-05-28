# Statistical model for denticle organization
#
# For more details, see [Spencer, Schaumberg and Zallen, Molecular Biology of the Cell, 2017](https://www.ncbi.nlm.nih.gov/pubmed/28404752). If this code has been helpful to you, please cite our paper.
#
# Basic methodology:
# The equation L = α*D + D(N-1) + α*D describes the general case for determining the distance between denticles (D) given the known cell length (L) and the number of denticles (N) in the cell, where α is the spacing ratio (the ratio of the average denticle-to-cell-edge distance to the average denticle-to-denticle distance).
#
# The ideal position of each denticle n_i can be described by the equation position_i = (α + n_i - 1) / (2α + N - 1) where α is the spacing ratio and N is the total number of denticles in the cell.
#
# This program generates simulated denticle positions for a dataset comprising (L,N,D) according to a series of truncated normal distributions with μ = the ideal position (calculated from the above equation) and σ = D/{a specified set of divisors}. Denticle spacing distances are calculated from the positional data. Simulated spacing distances are then compared to the in vivo spacing data.





from __future__ import division

from sys import argv
import numpy as np
import pandas as pd
import scipy.stats as sps
# import datetime as dt
import time
import string
import csv



# Functions

def IdealPositionMu(celllength, dentnumber, alpha):
    # calculates the ideal position(s) of denticle(s) in a cell with length = celllength
    distances = []
    for dent in range(1,dentnumber+1):
        distances.append(((alpha + dent - 1) / (2*alpha + dentnumber - 1)) * celllength)
    return distances


def AddByLine(datatoadd, addto_list, frame, idx):
    datatoadd = np.array(datatoadd)
    datatoadd.T
    for line in datatoadd:
        addto_list.append(tuple(line))

    tempseries = pd.Series(addto_list,index=idx)
    return frame.append(tempseries,ignore_index=True)


def SigmasCalculator(celllength, dentnumber, alphas, stdevnumbers):
    alphas = np.array(alphas)
    # stdevnumbers = np.array([[5], [6], [7], [8]])
    sigma_lengths = np.zeros((len(stdevnumbers),len(alphas)))
    sigma_lengths = (celllength / (2*alphas + dentnumber -1)) / stdevnumbers
    # sigmasFrame = pd.DataFrame(sigma_lengths)
    return sigma_lengths #, sigmasFrame


def GenerateDistributions(mus, sigmas, celllength):
    # generates values from a truncated normal distribution (adjusted to account for cell length - denticles cannot be placed outside of the cell boundaries, ie less than 0 or more than celllength)
    # mus is a tuple or list, sigmas is a series or dict, celllength and dentnumber are values

    positions = np.zeros((len(mus), len(sigmas)))

    for s, sigma in enumerate(sigmas):
        for m, mu in enumerate(mus):
            lower, upper = 0, celllength
            positions[m,s] = sps. ((lower-mu)/sigma, (upper-mu)/sigma, loc=mu, scale=sigma)

    positions = np.sort(positions,axis=0)
    # sort by rows to order denticle positions from smallest to largest
    distances = positions[1:,:] - positions[0:-1,:]

    return positions.T, distances.T

    # each row is an entire cell with all the denticles
    # each row of the returned arrays hold positions for a different value of mu, columns correspond to different sigmas

    # arr.T :                                             arr :
    #     [m1] [m2] [m3] [m4]                           #       [5] [6] [7] [8]
    # [5]   a   a   a   a                               # [m1]   a   a   a   a
    # [6]   b   b   b   b                               # [m2]   b   b   b   b
    # [7]   c   c   c   c                               # [m3]   c   c   c   c
    # [8]   d   d   d   d                               # [m4]   d   d   d   d

def SliceSaver(array,names):
    # takes a 3d array, slices by z axis, saves each slice as a csv
    for dex, arr in enumerate(array):
        np.savetxt(genotype +'_'+ names[dex] + '.csv', arr, delimiter=',')


def NameMaker(idx1, addon):
    names = []
    for i, elem in enumerate(idx1):
        names.append(str(elem) + '_' + addon + time.strftime("%Y%m%d_%H%M%S"))
        # for windows
        # names.append(str(elem) + '_' + addon)

    return names


def IndivStatTest(array, names, filename_out):
# IN: 3D np array, list of strings with length=arr[X,:,:] (array axis 0), name of csv file
    for dex, section in enumerate(array):
        for itnumber, iterationdata in enumerate(section.T):
            simdata = iterationdata[np.nonzero(iterationdata)]

            test_ks = sps.ks_2samp(invivo_d, simdata)
            test_mwu = sps.mannwhitneyu(invivo_d[np.nonzero(invivo_d)], simdata)

            with open(filename_out, 'a') as f:
                csv.writer(f).writerows([[itnumber, names[dex], test_ks[0], test_ks[1], TestPasses(test_ks[1], 0.05), test_mwu[0], test_mwu[1], TestPasses(test_mwu[1], 0.05)]])


def TestPasses(pval, cutoff):
    if pval <= cutoff:
        return 'different'
    elif pval > cutoff:
        return 'same'


def ComboStatTest(array,names,filename_out):
# IN: 3D np array, list of strings with length=arr[X,:,:] (array axis 0), name of csv file
    for dex, section in enumerate(array):
        simdata = section[np.nonzero(section)].flatten()

        test_ks = sps.ks_2samp(invivo_d, simdata)
        test_mwu = sps.mannwhitneyu(invivo_d, simdata)

        with open(filename_out, 'a') as f:
            csv.writer(f).writerows([['-', names[dex], test_ks[0], test_ks[1], TestPasses(test_ks[1], 0.05), test_mwu[0], test_mwu[1], TestPasses(test_mwu[1], 0.05)]])








starttime = time.clock()   # returns time in seconds

# for cmd line call
script, input_file, cells_to_simulate = argv
data = pd.read_csv(input_file)

# # for jupyter notebook
# input_file, cells_to_simulate = 'yw_CellbyCell.csv', 1
# data = pd.read_csv(input_file)


# get genotype from filename
genotype = input_file.split('_')[0]  # split the string at the underscore, returns a list of the parts (no underscores), take the first/before part


# clear out missing data
data = data.replace(0,np.nan)       # turn zeros into NaNs
data = data.dropna(how='all')       # drop any column (axis=0) or row (axis=1) where ALL values are NaN
data = data.replace(np.nan,0)       # turn NaNs back into zeros, so arithmetic can be normal for 'true-zero' values
data = data[data.dentincell != 1]   # drop columns where the cell has only 1 denticle

# save the invivo denticle separation distances to a 1xN matrix
invivo_d = data.as_matrix(columns=data.columns[10:]).flatten()
invivo_d = invivo_d[np.nonzero(invivo_d)]


# # using the 'absolute' DV length (dist between edge markers):
# invivo_ln = pd.DataFrame(data,columns=['dentincell','Dvlen'])
# invivo_ln.columns = ['dentincell','dvlength']

# using the summed DVlength (sum of dent-edge, dent-dent ... dent-dent, dent-edge); sum to get the additive, rather than absolute, DV length:
invivo_ln = pd.DataFrame(data,columns=['dentincell'])
invivo_ln['dvlength'] = data['dentEdgeL'] + data['dentEdgeR'] + (data[data.columns[10:]]).sum(axis=1,skipna=True,numeric_only=True)


# invivo_ln['dentincell'] = invivo_ln['dentincell'].astype(int)
invivo_ln = invivo_ln.T

# data = None       # clear data from memory




# Constants and basic indicies
maxDent = int(max(invivo_ln.ix['dentincell']))

alphas = (5/10, 6/10, 2/3, 7/10, 3/4, 8/10, 9/10, 1)
idx_alphas = ['halfD','sixtenthsD','twothirdsD','sevententhsD','threequartersD','eighttenthsD','ninetenthsD','oneD']

basicinfo = ['dvlength','dentincell']

idx_positions = basicinfo + idx_alphas


# define arrays for iteration data grouped by sigma value
#                   alphas           iterations         cells in dataset * max denticles
# sd_five     = np.zeros((len(idx_alphas), len(invivo_ln.T)*(int(maxDent)-1), int(cells_to_simulate)))
# sd_six      = np.zeros((len(idx_alphas), len(invivo_ln.T)*(int(maxDent)-1), int(cells_to_simulate)))
# sd_seven    = np.zeros((len(idx_alphas), len(invivo_ln.T)*(int(maxDent)-1), int(cells_to_simulate)))
# sd_eight    = np.zeros((len(idx_alphas), len(invivo_ln.T)*(int(maxDent)-1), int(cells_to_simulate)))


# numberofstdevs = np.array([[5], [6], [7], [8])
# idx_sigmas = ['five','six','seven','eight']
# # define arrays for iteration data grouped by sigma value
# #                   alphas           iterations         cells in dataset * max denticles
# sd_five     = np.zeros((len(idx_alphas), len(invivo_ln.T)*(int(maxDent)-1), int(cells_to_simulate)))
# sd_six      = np.zeros((len(idx_alphas), len(invivo_ln.T)*(int(maxDent)-1), int(cells_to_simulate)))
# sd_seven    = np.zeros((len(idx_alphas), len(invivo_ln.T)*(int(maxDent)-1), int(cells_to_simulate)))
# sd_eight    = np.zeros((len(idx_alphas), len(invivo_ln.T)*(int(maxDent)-1), int(cells_to_simulate)))
#
# sdarrs = [sd_five, sd_six, sd_seven, sd_eight]


numberofstdevs = np.array([[5], [6], [7], [8], [9]])
idx_sigmas = ['five','six','seven','eight','nine']
# sd_five     = np.zeros((len(idx_alphas), len(invivo_ln.T)*(int(maxDent)-1), int(cells_to_simulate)))
# sd_six      = np.zeros((len(idx_alphas), len(invivo_ln.T)*(int(maxDent)-1), int(cells_to_simulate)))
# sd_seven    = np.zeros((len(idx_alphas), len(invivo_ln.T)*(int(maxDent)-1), int(cells_to_simulate)))
# sd_eight    = np.zeros((len(idx_alphas), len(invivo_ln.T)*(int(maxDent)-1), int(cells_to_simulate)))
# sd_nine     = np.zeros((len(idx_alphas), len(invivo_ln.T)*(int(maxDent)-1), int(cells_to_simulate)))

# sdarrs = [sd_five, sd_six, sd_seven, sd_eight, sd_nine]




# MAIN

# initial calculations of the ideal positions and sigmas
positionframe = pd.DataFrame()
sigmaframe = pd.DataFrame()

for dex, row in enumerate(invivo_ln):
    celllength = invivo_ln[row]['dvlength']
    dentnumber = invivo_ln[row]['dentincell'].astype(int)

    # calculate mu (ideal positions); generates a len(invivo_ln) x 10 DataFrame
    positionlist = []
    for alpha in alphas:
        positionlist.append(IdealPositionMu(celllength,dentnumber,alpha))
    positionframe = AddByLine(positionlist, [celllength, dentnumber], positionframe, idx_positions)

    # calculate sigmas; generates a multi-indexed len(invivo_ln)*len(sigmas) x 8 DataFrame
    sigmaframe = pd.concat([sigmaframe,
                                pd.DataFrame(SigmasCalculator(celllength,dentnumber,alphas,numberofstdevs),
                                index=[[dex]*len(idx_sigmas),idx_sigmas],columns=idx_alphas)])

positionframe = positionframe[idx_positions]

sigmaframe.index.names = ['cell','numbSDs']
sigmaframe.columns.names = ['alpha']
positionframe.index.names = ['cell']

# positionframe.to_csv('ideal_positions_' + dt.datetime.isoformat(dt.datetime.now()) + '_' + str(time.clock() - starttime) + '.csv')
# sigmaframe.to_csv('sigmavalues_' + dt.datetime.isoformat(dt.datetime.now()) + '_' + str(time.clock() - starttime) + '.csv')

positionframe.to_csv('ideal_positions.csv')
sigmaframe.to_csv('sigmavalues.csv')


# print('initial calcs',time.clock() - starttime)

initcalctime = time.clock() - starttime





# generate random positions


for iteration in range(0, int(cells_to_simulate)):
    ittime = time.clock()
    #                           number of cells    each (alpha,sigma) pair + L + N         denticles possible
    randompositions = np.zeros((len(invivo_ln.T), (len(idx_sigmas) * len(idx_alphas) + 2), maxDent))
    randomdistances = np.zeros((len(invivo_ln.T), (len(idx_sigmas) * len(idx_alphas) + 2), maxDent))

    for cellnumber, __ in enumerate(positionframe.T):
        celllength, dentnumber = positionframe.iloc[cellnumber]['dvlength'], int(positionframe.iloc[cellnumber]['dentincell'])

        randomdistances[cellnumber, 0, :] = positionframe.iloc[cellnumber]['dvlength']
        randomdistances[cellnumber, 1, :] = positionframe.iloc[cellnumber]['dentincell']

        randompositions[cellnumber, 0, :] = positionframe.iloc[cellnumber]['dvlength']
        randompositions[cellnumber, 1, :] = positionframe.iloc[cellnumber]['dentincell']


        for alphano, elem in enumerate(idx_alphas):
            positionset = positionframe.iloc[cellnumber][elem]
            sigmaset = sigmaframe.ix[positionframe.iloc[cellnumber].name][elem]

            positions, distances = GenerateDistributions(positionframe.iloc[cellnumber][elem], sigmaframe.ix[positionframe.iloc[cellnumber].name][elem], celllength)
            randompositions[cellnumber, (len(idx_sigmas)*alphano)+len(basicinfo) : (len(idx_sigmas)*alphano + len(idx_sigmas))+len(basicinfo), 0:dentnumber] = positions
            randomdistances[cellnumber, (len(idx_sigmas)*alphano)+len(basicinfo) : (len(idx_sigmas)*alphano + len(idx_sigmas))+len(basicinfo), 0:dentnumber-1] = distances

            # sd_five[alphano, (dentnumber-1)*cellnumber : (dentnumber-1)*cellnumber+(dentnumber-1), iteration] = distances[0,:]
            # sd_six[alphano, (dentnumber-1)*cellnumber : (dentnumber-1)*cellnumber+(dentnumber-1), iteration] = distances[1,:]
            # sd_seven[alphano, (dentnumber-1)*cellnumber : (dentnumber-1)*cellnumber+(dentnumber-1), iteration] = distances[2,:]
            # sd_eight[alphano, (dentnumber-1)*cellnumber : (dentnumber-1)*cellnumber+(dentnumber-1), iteration] = distances[3,:]


    # print('iteration time',time.clock() - ittime)

    # items: axis 0, each item corresponds to a DataFrame contained inside; major_axis: axis 1, it is the index (rows) of each of the DataFrames; minor_axis: axis 2, it is the columns of each of the DataFrames
    randompositionpanel = pd.Panel(randompositions,
                            major_axis=[list(np.repeat('basicinfo',2)) + list(np.repeat(idx_alphas, len(idx_sigmas))), ['celllength','dentnumber'] + list(idx_sigmas * len(idx_alphas))],
                            minor_axis=list(string.ascii_lowercase[:maxDent]))

    randomdistancepanel = pd.Panel(randomdistances,
                            major_axis=[list(np.repeat('basicinfo',2)) + list(np.repeat(idx_alphas, len(idx_sigmas))), ['celllength','dentnumber'] + list(idx_sigmas * len(idx_alphas))],
                            minor_axis=list(string.ascii_lowercase[:maxDent]))

    randompositionframe = randompositionpanel.to_frame()
    randomdistanceframe = randomdistancepanel.to_frame()


    elapsedtime = time.clock() - ittime

    randompositionframe.to_csv(genotype + '_iteration' + str(iteration) + '_random_positions_' + time.strftime("%Y%m%d_%H%M%S") + '_iterationtime=' + str(elapsedtime) + '.csv')
    randomdistanceframe.to_csv(genotype + '_iteration' + str(iteration) + '_random_distances_' + time.strftime("%Y%m%d_%H%M%S") + '_iterationtime=' + str(elapsedtime) + '.csv')



# for dex, varname in enumerate(sdarrs):

#     SliceSaver(varname,NameMaker(idx_alphas,idx_sigmas[dex]))

#     IndivStatTest(varname, NameMaker(idx_alphas, idx_sigmas[dex]), genotype + '_statistics_' + str(dex) + idx_sigmas[dex] + time.strftime("%Y%m%d_%H%M%S") + '.csv')
#     ComboStatTest(varname, NameMaker(idx_alphas, idx_sigmas[dex]), genotype + '_statistics_' + str(dex) + idx_sigmas[dex] + time.strftime("%Y%m%d_%H%M%S") + '.csv')

    # # for windows
    # IndivStatTest(varname, NameMaker(idx_alphas, idx_sigmas[dex]), genotype + '_statistics_' + str(dex) + idx_sigmas[dex] + '.csv')
    # ComboStatTest(varname, NameMaker(idx_alphas, idx_sigmas[dex]), genotype + '_statistics_' + str(dex) + idx_sigmas[dex] + '.csv')


endtime = time.clock()



# # save data from each alpha,sigma combo as a separate file
# SliceSaver(sd_five,NameMaker(idx_alphas,'five'))
# SliceSaver(sd_six,NameMaker(idx_alphas,'six'))
# SliceSaver(sd_seven,NameMaker(idx_alphas,'seven'))
# SliceSaver(sd_eight,NameMaker(idx_alphas,'eight'))


# IndivStatTest(sd_five, NameMaker(idx_alphas, 'five'), genotype + '_statistics_five.csv')
# ComboStatTest(sd_five, NameMaker(idx_alphas, 'five'), genotype + '_statistics_five.csv')

# IndivStatTest(sd_six, NameMaker(idx_alphas, 'six'), genotype + '_statistics_six.csv')
# ComboStatTest(sd_six, NameMaker(idx_alphas, 'six'), genotype + '_statistics_six.csv')

# IndivStatTest(sd_seven, NameMaker(idx_alphas, 'seven'), genotype + '_statistics_seven.csv')
# ComboStatTest(sd_seven, NameMaker(idx_alphas, 'seven'), genotype + '_statistics_seven.csv')

# IndivStatTest(sd_eight, NameMaker(idx_alphas, 'eight'), genotype + '_statistics_eight.csv')
# ComboStatTest(sd_eight, NameMaker(idx_alphas, 'eight'), genotype + '_statistics_eight.csv')




# print('iteration time',elapsedtime)
# print('total time',endtime-starttime)

