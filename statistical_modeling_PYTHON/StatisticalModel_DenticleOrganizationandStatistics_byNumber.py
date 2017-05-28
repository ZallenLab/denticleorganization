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
import glob
import os
import itertools
import time
import string
import csv

import numpy as np
import pandas as pd
import scipy.stats as sps

import seaborn
import matplotlib.pyplot as plt



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
    # each row is an entire cell with all the denticles
    # each row of the returned arrays hold positions for a different value of mu, columns correspond to different sigmas
    # arr.T :                      # arr :
    #     [m_1] [m_2] [m_3]        #       [5] [6] [7]
    # [5]   a     a     a          # [m1]   a   a   a
    # [6]   b     b     b          # [m2]   b   b   b
    # [7]   c     c     c          # [m3]   c   c   c

    positions = np.zeros((len(mus), len(sigmas)))

    for s, sigma in enumerate(sigmas):
        for m, mu in enumerate(mus):
            lower, upper = 0, celllength
            positions[m,s] = sps.truncnorm.rvs((lower-mu)/sigma, (upper-mu)/sigma, loc=mu, scale=sigma)

    positions = np.sort(positions,axis=0)
    # sort by rows to order denticle positions from smallest to largest
    distances = positions[1:,:] - positions[0:-1,:]
    # calculate the distance between the points

    return positions.T, distances.T


def IndivStatTest(simdata, model, sigma, filename_out):
# IN: 3D np array, list of strings with length=arr[X,:,:] (array axis 0), name of csv file
    simdata = simdata[np.nonzero(simdata)]

    test_ks = sps.ks_2samp(invivo_d, simdata)
    # outputs [ks-statistic, p-value]
    # If the K-S statistic is small or the p-value is high, then we cannot reject the hypothesis that the distributions of the two samples are the same.
    test_mwu = sps.mannwhitneyu(invivo_d[np.nonzero(invivo_d)], simdata)
    # returns [mwu-statistic, p-value];  One-sided p-value assuming a asymptotic normal distribution.
    # Use only when the number of observation in each sample is > 20 and you have 2 independent samples of ranks. Mann-Whitney U is significant if the u-obtained is LESS THAN or equal to the critical value of U.
    # This test corrects for ties and by default uses a continuity correction. The reported p-value is for a one-sided hypothesis, to get the two-sided p-value multiply the returned p-value by 2.

    with open(filename_out, 'a') as f:
        csv.writer(f).writerows([[model+'_'+sigma,  test_ks[0], test_ks[1], TestPasses(test_ks[1], 0.05), test_mwu[0], test_mwu[1], TestPasses(test_mwu[1], 0.05)]])


def TestPasses(pval, cutoff):
    if pval < cutoff:
        return 'different'
    elif pval >= cutoff:
        return 'same'











starttime = time.clock()   # returns time in seconds

# for cmd line call
script, input_file, cells_to_simulate = argv
dataIN = pd.read_csv(input_file)

basepath = os.getcwd()

# # for jupyter notebook
# input_file, cells_to_simulate = 'yw_CellbyCell.csv', 1
# dataIN = pd.read_csv(input_file)


# get genotype from filename
genotype = input_file.split('/')[-1].split('_')[-2]
# genotype = input_file.split('_')[-1]  # split the string at the underscore, returns a list of the parts (no underscores), take the first/before part
# genotype = 'Larvae'
print(genotype)

# clear out missing dataIN
dataIN = dataIN.replace(0,np.nan)       # turn zeros into NaNs
dataIN = dataIN.dropna(how='all')       # drop any column (axis=0) or row (axis=1) where ALL values are NaN
dataIN = dataIN.replace(np.nan,0)       # turn NaNs back into zeros, so arithmetic can be normal for 'true-zero' values
# data = data[data.dentincell != 1]   # drop columns where the cell has only 1 denticle
for dn in range(2, int(max(dataIN.dentincell.unique())) + 1):

    data = dataIN[dataIN.dentincell == dn]   # drop columns where the cell has only 1 denticle
    numberforfilenames = '_' + str(dn) + 'denticlesonly_'

    # save the invivo denticle separation distances to a 1xN matrix
    invivo_d = data.as_matrix(columns=data.columns[10:]).flatten()
    invivo_d = invivo_d[np.nonzero(invivo_d)]


    # using the 'absolute' DV length (dist between edge markers):
    invivo_ln = pd.DataFrame(data,columns=['dentincell','Dvlen'])
    invivo_ln.columns = ['dentincell','dvlength']
    numberforfilenames = numberforfilenames + 'AbsoluteDV_'

    # # using the summed DVlength (sum of dent-edge, dent-dent ... dent-dent, dent-edge); sum to get the additive, rather than absolute, DV length:
    # invivo_ln = pd.DataFrame(data,columns=['dentincell'])
    # invivo_ln['dvlength'] = data['dentEdgeL'] + data['dentEdgeR'] + (data[data.columns[10:]]).sum(axis=1,skipna=True,numeric_only=True)
    # numberforfilenames = numberforfilenames + 'SummedDV_'

    invivo_ln = invivo_ln.T






    # Constants and basic indicies
    maxDent = int(max(invivo_ln.ix['dentincell']))

    alphas = (5/10, 6/10, 2/3, 7/10, 3/4, 8/10, 9/10, 1)
    idx_alphas = ['halfD','sixtenthsD','twothirdsD','sevententhsD','threequartersD','eighttenthsD','ninetenthsD','oneD']
    basicinfo = ['dvlength','dentincell']
    idx_positions = basicinfo + idx_alphas

    numberofstdevs = np.array([[5], [6], [7], [8], [9]])
    idx_sigmas = ['five','six','seven','eight','nine']






    # MAIN

    # Now that fundamentals are in place (data needed is imported from the base dir) make a new dir for this run, and move to it, so all the output files get placed here not the main dir
    dirname = genotype +'_Output' + numberforfilenames + time.strftime("%Y%m%d_%H%M%S")

    os.mkdir(dirname)
    os.chdir(dirname)




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


    positionframe.to_csv(genotype + numberforfilenames + 'ideal_positions' + '.csv')
    sigmaframe.to_csv(genotype + numberforfilenames + 'sigmavalues' + '.csv')


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

        randompositionframe.to_csv(genotype + '_iteration' + str(iteration) + '_random_positions' + numberforfilenames + time.strftime("%Y%m%d_%H%M%S") + '_iterationtime=' + str(elapsedtime) + '.csv')
        randomdistanceframe.to_csv(genotype + '_iteration' + str(iteration) + '_random_distances' + numberforfilenames + time.strftime("%Y%m%d_%H%M%S") + '_iterationtime=' + str(elapsedtime) + '.csv')


    endtime = time.clock()
    # end of simulation-data-generation steps











    # calculate statistics of datasets

    # establish basics again
    # TODO this can be merged with the initial things above 'MAIN'
    # idx_alphas =  ['halfD','sixtenthsD','twothirdsD','sevententhsD','threequartersD','eighttenthsD','ninetenthsD','oneD']
    # idx_sigmas = ['five','six','seven','eight','nine']


    #
    # example filename: yw_iteration9940_random_distances_20160127_192532_iterationtime=24.116169999993872.csv
    names = []
    for fname in glob.glob('*iteration*distances*.csv'):
        names.append(fname)

    now = time.strftime("%Y%m%d_%H%M%S")


    for name in names:
        frame = pd.read_csv(name, index_col=(0,1,2))
        # read in csv file, with columns indices specificed
        frame = frame.T
        frame.columns.names = [None, None, None]
        # flip, because it's easier to work with columns than with rows indexed... and clear out names

        for elem in itertools.product(idx_alphas,idx_sigmas):
            model,sigma = elem
            subframe = frame[model][sigma].values.flatten()
            csvname = genotype + numberforfilenames + sigma + '_statistics_' + now +'.csv'
            IndivStatTest(subframe, model, sigma, csvname)





    # summarize statistics data
    names = []
    for fname in glob.glob('*_statistics_*.csv'):
        names.append(fname)

    idx_stats = ['iter','model','ks_stat','ks_pvalue','ks_sameness','mwu_stat','mwu_pvalue','mwu_sameness']
    terms = {'same':1,'different':0}
    splitter = 'ks'

    for name in names:
        genotype = name.split('_')[0]
        frame = pd.read_csv(name,header=None, names=idx_stats[1:])
        frame = frame.sort_values(by='model')

        frame.index = list(range(0,int(len(frame)/len(idx_alphas))))*len(idx_alphas)

        stdev = name.split('_')[3]

        fr2 = pd.DataFrame()
        fr3 = pd.DataFrame()
        fr4 = pd.DataFrame()

        for modelname in (model + '_' + stdev for model in idx_alphas):
            temp = frame[frame['model'] == modelname]

            fr2[modelname.split('2')[0].split('_')[0]] = frame[frame['model'] == modelname][splitter + '_pvalue']
            fr3[modelname.split('2')[0].split('_')[0]] = frame[frame['model'] == modelname][splitter + '_sameness']
            fr3[modelname.split('2')[0].split('_')[0]] = fr3[modelname.split('2')[0].split('_')[0]].map(terms)

        # save
        fr2.to_csv(genotype + numberforfilenames + stdev + '_' + splitter + '_pvalue' + '.csv')
        fr3.to_csv(genotype + numberforfilenames + stdev + '_' + splitter + '_sameness' + '.csv')

        # calculate fractions that pass/fail, summarize, and make plots
        fr4['p<0.05'] = (1-(fr3.sum()/len(fr3)))*100
        fr4['p>0.05'] = (fr3.sum()/len(fr3))*100

        myplot = fr4.plot(kind='bar',stacked=True,color=['b','y'])
        myplot.legend(bbox_to_anchor=(.55,1.1),loc=2,ncol=2,borderaxespad=0.5)
        plt.title('a = ' + modelname.split('2')[0].split('_')[0] + ', s = ' + stdev, loc='left')
        fig = myplot.get_figure()
        fig.savefig(genotype + numberforfilenames + stdev + '.png', bbox_inches='tight')
        plt.close(fig)

        fr4.to_csv(genotype + numberforfilenames + stdev + '_' + splitter + '_fractions' + '.csv')


    os.chdir(basepath)

