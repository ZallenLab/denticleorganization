
# coding: utf-8

# In[1]:

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


# In[2]:
def IndivStatTest(original_data,test_data, model, sigma, filename_out):
# IN: 3D np array, list of strings with length=arr[X,:,:] (array axis 0), name of csv file
    test_data = test_data[np.nonzero(test_data)]

    test_ks = sps.ks_2samp(original_data, test_data)
    # outputs [ks-statistic, p-value]
    # If the K-S statistic is small or the p-value is high, then we cannot reject the hypothesis that the distributions of the two samples are the same.
    test_mwu = sps.mannwhitneyu(original_data[np.nonzero(original_data)], test_data)
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




#  Import 

starttime = time.clock()   # returns time in seconds

# # for cmd line call
script, genotype, simnumber_to_test, val_alpha, val_sigma = argv

os.chdir(glob.glob('yw*/')[0])

input_file = glob.glob('*iteration' + simnumber_to_test + '_random_distances*.csv')[0]
simulation_data = pd.read_csv(input_file, index_col=(0, 1, 2))


# # for jupyter notebook
# simnumber_to_test = '1'
# input_file = glob.glob('*iteration' + simnumber_to_test + '_random_distances*.csv')[0]
# simulation_data = pd.read_csv(input_file, index_col=(0,1,2))

# val_alpha = 'twothirdsD'
# val_sigma = 'six'


# create saving path 
dirname = 'Validation_against_'+ val_alpha +'_'+ val_sigma +'_iternumber_' + simnumber_to_test +'_'+ genotype + time.strftime("%Y%m%d_%H%M%S")
os.mkdir(dirname)

savepath = os.path.join(os.getcwd(),dirname)


# get genotype from filename
# genotype = input_file.split('_')[0]  # split the string at the underscore, returns a list of the parts (no underscores), take the first/before part

simuldata = simulation_data.T
simuldata = simuldata[val_alpha][val_sigma].values.flatten()
simuldata = simuldata[np.nonzero(simuldata)]


# Constants and basic indicies

alphas = (5/10, 6/10, 2/3, 7/10, 3/4, 8/10, 9/10, 1)
idx_alphas = ['halfD', 'sixtenthsD', 'twothirdsD', 'sevententhsD', 'threequartersD', 'eighttenthsD', 'ninetenthsD', 'oneD']
basicinfo = ['dvlength','dentincell']
idx_positions = basicinfo + idx_alphas

numberofstdevs = np.array([[5], [6], [7], [8], [9]])
idx_sigmas = ['five', 'six', 'seven', 'eight', 'nine']


# calculate statistics of datasets

# establish basics again
# TODO this can be merged with the initial things above 'MAIN'
# idx_alphas =  ['halfD','sixtenthsD','twothirdsD','sevententhsD','threequartersD','eighttenthsD','ninetenthsD','oneD']
# idx_sigmas = ['five','six','seven','eight','nine']


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
        csvname = 'Validation_against_'+ val_alpha +'_'+ val_sigma +'_' + genotype +'_statistics_' + sigma +'_'+ now +'.csv'
        pth = os.path.join(savepath,csvname)
        IndivStatTest(simuldata, subframe, model, sigma, pth)





# summarize statistics data 
os.chdir(savepath)

names = []
for fname in glob.glob('Validation_against_*_statistics_*.csv'):
    names.append(fname)

idx_stats = ['iter','model','ks_stat','ks_pvalue','ks_sameness','mwu_stat','mwu_pvalue','mwu_sameness']
terms = {'same':1,'different':0}
splitter = 'ks'

for name in names: 
#     genotype = name.split('_')[0]
    frame = pd.read_csv(name,header=None, names=idx_stats[1:])
    frame = frame.sort_values(by='model')

    frame.index = list(range(0,int(len(frame)/len(idx_alphas))))*len(idx_alphas)
    
    stdev = name.split('_')[6]

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
    fr2.to_csv(os.path.join(savepath, genotype + '_' + stdev + '_' + splitter + '_pvalue' + '.csv'))
    fr3.to_csv(os.path.join(savepath, genotype + '_' + stdev + '_' + splitter + '_sameness' + '.csv'))

    # calculate fractions that pass/fail, summarize, and make plots
    fr4['p<0.05'] = (1-(fr3.sum()/len(fr3)))*100
    fr4['p>0.05'] = (fr3.sum()/len(fr3))*100

    myplot = fr4.plot(kind='bar',stacked=True,color=['b','y'])
    myplot.legend(bbox_to_anchor=(.55,1.1),loc=2,ncol=2,borderaxespad=0.5)
    plt.title('s = ' + stdev, loc='left')
    fig = myplot.get_figure()
    fig.savefig(os.path.join(savepath, 'Validation_against_'+ val_alpha +'_'+ val_sigma +'_' + genotype+'_'+stdev+'.png'),bbox_inches='tight')
#     plt.close(fig)

    fr4.to_csv(os.path.join(savepath, 'Validation_against_'+ val_alpha +'_'+ val_sigma +'_' + genotype + '_' + stdev + '_' + splitter + '_fractions' + '.csv'))

