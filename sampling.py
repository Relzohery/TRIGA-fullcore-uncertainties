from scipy import interpolate
import numpy as np
import os
import json
from compositions import get_compositions
from Newt_format import NewtInput, Isotopes_of_interest


Uncertainty_data = {}
config = 'Config.txt'  # text file of the core configuration
locs = {}
with open (config, 'r') as f:
    lines = f.readlines()
    for line in lines:
        k = line.split()
        # check the fuel locations, and store them
        if k[1].isdigit():
            locs[k[1]] = k[0]


UQ_data ={}
Bu_data = 'BU_uncertainty.txt' #text file contains fuel elements, burnups and uncertainties
with open(Bu_data) as f:
     lines = f.readlines()[1:]
     for line in lines:
        k = line.split()
        # check if the fuel element in the current core configuration.
        if k[0] in locs.keys():
            UQ_data[k[0]] = (float(k[1]), float(k[2]))



if __name__ == '__main__':
    # produce the reference case without perturbations
    core_comp = {}
    for element, data in UQ_data.items():
        BU = float(data[0])
        core_comp[element] = get_compositions(BU).comp
    comp_file = './core_comp.json'
    with open(comp_file, 'w') as f:
        json.dump(core_comp, f)
    Input = NewtInput('Config.txt', Isotopes_of_interest, core_comp = comp_file,
                      Inputname = 'Triga_2D_Unperturbed.inp')
    Input.make_input()

    # produce 100 perturbed samples to be used for forward sampling UQ
    no_samples = 100
    for sample in range(no_samples):
        try:
            # make folder for each sample
            os.mkdir('sample{}'.format(sample))
        except:
            pass
        core_comp = {}
        for element, data in UQ_data.items():
        #for loc, data in covariance.items():
            #sampling the burnup
            mean = data[0]
            # uncertainties are given as a percentage so they have to be
            # converted to their absolute values
            std = data[0]*data[1]/100
            BU = np.random.normal(mean, std)
            core_comp[element]= get_compositions(BU).comp
        comp_file = './sample{}/core_comp.json'.format(sample)
        with open(comp_file, 'w') as f:
            json.dump(core_comp, f)
        Input = NewtInput('Config.txt', Isotopes_of_interest, core_comp = comp_file,
                          Inputname = './sample{}/sample{}.inp'.format(sample, sample))
        # write NewtInput
        Input.make_input()

