# required command-line arguments: <in vivo dataset file>  <number of iterations>
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
import num2words





# Functions
def GenerateDist_Uniformrandom(celllength, denticlenumber):
    positions = np.random.rand(denticlenumber)
    positions = np.sort(positions,axis=0)
    # sort smallest to largest
    relativedistances = positions[1:] - positions[0:-1]
    absolutedistances = relativedistances * celllength
    # calculate the distance between the points

    return positions, relativedistances, absolutedistances

def IndivStatTest(simdata, model, filename_out):
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
        csv.writer(f).writerows([[model,  test_ks[0], test_ks[1], TestPasses(test_ks[1], 0.05), test_mwu[0], test_mwu[1], TestPasses(test_mwu[1], 0.05)]])


def TestPasses(pval, cutoff):
    if pval < cutoff:
        return 'different'
    elif pval >= cutoff:
        return 'same'









starttime = time.clock()   # returns time in seconds

# for cmd line call
script, input_file, cells_to_simulate = argv
data = pd.read_csv(input_file)

# for jupyter notebook
# input_file, cells_to_simulate = 'yw_CellbyCell.csv', 1
# data = pd.read_csv(input_file)


# get genotype from filename
genotype = input_file.split('_')[0]     # split the string at the underscore, returns a list of the parts (no underscores), take the first/before part


# clear out missing data
data = data.replace(0, np.nan)       # turn zeros into NaNs
data = data.dropna(how='all')   # drop any column (axis=0) or row (axis=1) where ALL values are NaN
data = data.replace(np.nan, 0)  # turn NaNs back into zeros, so arithmetic can be normal for 'true-zero' values
data = data[data.dentincell != 1]   # drop columns where the cell has only 1 denticle

# save the invivo denticle separation distances to a 1xN matrix
invivo_d = data.as_matrix(columns=data.columns[10:]).flatten()
invivo_d = invivo_d[np.nonzero(invivo_d)]


# using the 'absolute' DV length (dist between edge markers):
dvlen_type = 'AbsoluteDVlength'
invivo_ln = pd.DataFrame(data,columns=['dentincell','Dvlen'])
invivo_ln.columns = ['dentincell','dvlength']

# # using the summed DVlength (sum of dent-edge, dent-dent ... dent-dent, dent-edge); sum to get the additive, rather than absolute, DV length:
# dvlen_type = 'SummedDVlength'
# invivo_ln = pd.DataFrame(data, columns=['dentincell'])
# invivo_ln['dvlength'] = data['dentEdgeL'] + data['dentEdgeR'] + (data[data.columns[10:]]).sum(axis=1, skipna=True, numeric_only=True)

invivo_ln = invivo_ln.T





# set dir and make new folder
dirname = genotype + '_Output_MonteCarlo_' + dvlen_type + '_' + time.strftime("%Y%m%d_%H%M%S")
os.mkdir(dirname)
os.chdir(dirname)




maxDent = int(max(invivo_ln.ix['dentincell']))


# generate random positions
# maxDent = int(max(invivo_ln.ix['dentincell']))
# maxDent = int(max(invivo_ln.ix['dentincell']))

# randompositions = np.zeros((len(invivo_ln.T), (cells_to_simulate), maxDent))
# randomdistances_rel = np.zeros((len(invivo_ln.T), (cells_to_simulate), (maxDent-1)))
# randomdistances_abs = np.zeros((len(invivo_ln.T), (cells_to_simulate), (maxDent-1)))

# randompositions[:,0:2,0] = invivo_ln.T

info = invivo_ln.T


# generate random positions
for iteration in range(0, int(cells_to_simulate)):

    ittime = time.clock()

    positions = np.zeros((len(invivo_ln.T),maxDent))

    for dex, cellinfo in enumerate(invivo_ln.T.values):
        denticlenumber, celllength = cellinfo
        denticlenumber = int(denticlenumber)

        positions[dex,:denticlenumber] = np.sort(np.random.rand(int(denticlenumber)),axis=0)

    # sort smallest to largest
    relativedistances = positions[:,1:] - positions[:, 0:-1]
    relativedistances[relativedistances<0] = 0
    # calculate the distance between the points
    absolutedistances =(info.dvlength.values * relativedistances.T).T

    IndivStatTest(absolutedistances.reshape(absolutedistances.size,1), 'AbsoluteDistances','MonteCarlo_absdist_statistics.csv')

    adframe = pd.DataFrame(np.concatenate([info.values, absolutedistances], axis=1),columns=[['dentincell','celllength']+list(string.ascii_lowercase[:maxDent-1])])
    adframe['dentword'] = adframe.dentincell.map(lambda x: num2words.num2words(x))
    adframe = adframe.drop('dentincell',axis=1)
    adframe = adframe.replace(0,np.nan)
    adframe = adframe.set_index(['dentword','celllength'])
    adframe = adframe.sort_index(axis=0)
    adframe = adframe.stack()



    elapsedtime = time.clock() - ittime

    np.savetxt(genotype +'_'+str(iteration) +'_relativepositions_' + dvlen_type +'_'+ time.strftime("%Y%m%d_%H%M%S") + '_iterationtime=' + str(elapsedtime) + '.csv',positions,delimiter=',')
    np.savetxt(genotype +'_'+str(iteration) +'_'+'_relativedistances_' + dvlen_type +'_'+ time.strftime("%Y%m%d_%H%M%S") + '_iterationtime=' + str(elapsedtime) + '.csv',relativedistances,delimiter=',')



    adframe.to_csv(genotype +'_'+ str(iteration) + '_MonteCarlo_absolute_distances_' + dvlen_type +'_'+ time.strftime("%Y%m%d_%H%M%S") + '_iterationtime=' + str(elapsedtime) + '.csv')



    for number in (num2words.num2words(x) for x in range(1,maxDent)):
        dentnumb = adframe[number]
        dentnumb.to_csv(genotype+'_'+number+'_'+str(iteration) + '_absolute_distances.csv')











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
        csvname = genotype +'_statistics_' + sigma +'_'+ now +'.csv'
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

    stdev = name.split('_')[2]

    fr2 = pd.DataFrame()
    fr3 = pd.DataFrame()
    fr4 = pd.DataFrame()

#     for modelname in (y+z for y,z in itertools.product((x+'_' for x in idx_alphas),stdev)):
    for modelname in (model + '_' + stdev for model in idx_alphas):
        temp = frame[frame['model'] == modelname]

        fr2[modelname.split('2')[0].split('_')[0]] = frame[frame['model'] == modelname][splitter + '_pvalue']
        fr3[modelname.split('2')[0].split('_')[0]] = frame[frame['model'] == modelname][splitter + '_sameness']
        fr3[modelname.split('2')[0].split('_')[0]] = fr3[modelname.split('2')[0].split('_')[0]].map(terms)

    # save
    fr2.to_csv(genotype + '_' + stdev + '_' + splitter + '_pvalue' + '.csv')
    fr3.to_csv(genotype + '_' + stdev + '_' + splitter + '_sameness' + '.csv')

    # calculate fractions that pass/fail, summarize, and make plots
    fr4['p<0.05'] = (1-(fr3.sum()/len(fr3)))*100
    fr4['p>0.05'] = (fr3.sum()/len(fr3))*100

    myplot = fr4.plot(kind='bar',stacked=True,color=['b','y'])
    myplot.legend(bbox_to_anchor=(.55,1.1),loc=2,ncol=2,borderaxespad=0.5)
    plt.title('a = ' + modelname.split('2')[0].split('_')[0] + ', s = ' + stdev, loc='left')
    fig = myplot.get_figure()
    fig.savefig(genotype+'_'+stdev+'.png',bbox_inches='tight')
    plt.close(fig)

    fr4.to_csv(genotype + '_' + stdev + '_' + splitter + '_fractions' + '.csv')



