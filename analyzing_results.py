# -*- coding: utf-8 -*-
"""
Created on Sun Mar  4 20:40:24 2018

@author: R A B A B
"""
import numpy as np
from sampling import locs
import scipy.io as spio
import matplotlib.pylab as plt
import pickle

no_snapshots = 100
with open('fullcore_newt_results.p','rb') as f:
    data =  pickle.load(f)
iso = ['U235', 'U238', 'Pu239', 'Am241']

pop_mean = data['unperturbed']
results = data['results']
samples = data['samples']

# statistics analysis

sample_mean = np.mean(results, axis = 3)
std = np.std(results, axis = 3)


locations = sorted(locs.values())
#### plotting ###
#pop_means vs mean
plt.plot(range(no_snapshots), sample_mean[67, 10, :])
plt.plot(range(no_snapshots), pop_mean[67, 10, :])


# plotting uncertainites
plt.figure(2)
plt.errorbar(range(no_snapshots), sample_mean[4,3, :], std[4, 3, :])


###input _samples""