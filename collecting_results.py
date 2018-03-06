#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  4 10:45:26 2018

@author: rabab
"""

import numpy as np
from Newt_format import  Isotopes_of_interest
from sampling import locs, UQ_data
import scipy.io as spio
import json
import pickle

no_samples = 100
no_snapshots = 100

########## No perturbation case ###########

pop_mean = np.zeros((len(locs),len(Isotopes_of_interest), no_snapshots ))
fname = './RefCase/Triga_2D_Unperturbed.0000000000000000{}.plt'
for loc in range(1, len(locs)+1):
        name = '0' + str(loc) if loc<10 else loc
        pop_mean[loc-1, :, :] = np.loadtxt(fname.format(name), skiprows = 6, usecols = range(101)[1:])[:-2]


fname = './samples/sample{}/sample{}.0000000000000000{}.plt'
results = np.zeros((len(locs), len(Isotopes_of_interest), no_snapshots, no_samples))
samples = np.zeros(results.shape[:3]) # the samples are the isotopes at time=0

# collecting the samples for each core location
for sample in range(no_samples):
    for loc in range(1, len(locs)+1):
        name = '0' + str(loc) if loc<10 else loc
        results[loc-1, :, :, sample] = np.loadtxt(fname.format(sample,sample, name), skiprows = 6, usecols = range(101)[1:])[:-2]
        samples[loc-1, :, sample] = results[loc-1, :, 0,  sample]


fname = './samples/sample{}/core_comp.json'
#
samples_in = {}
for element in locs.keys():
    samples_in[locs[element]] = {}
    for Iso in Isotopes_of_interest:
        samples_in[locs[element]][Iso] = np.zeros(no_samples)

for sample in range(no_samples):
    with open(fname.format(sample)) as f:
        comp = json.load(f)
        for element in comp.keys():
            for Iso in Isotopes_of_interest:
                samples_in[locs[element]][Iso][sample]= comp[element][Iso]


with open ('fullcore_newt_results.p','wb') as f:
    pickle.dump({'results': results, 'unperturbed':pop_mean,'samples': samples, 'elements_data' :UQ_data }, f)



