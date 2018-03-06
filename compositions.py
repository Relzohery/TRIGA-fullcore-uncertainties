import pickle
from scipy import interpolate

# File containing the very recent data , 61 isotopes, 3 axial segments, 151 time step.
comp_loc = 'serp_total_BU_pic_600.p'
# File containing old data produced by Saqr, 31 isotopes, 7 axial segments, 41 time step.
#comp_loc = '1.p'


class get_compositions():
    ''' this calss takes a burnup value ans return the  a dictionary contains the corresponding compositions'''

    def __init__(self, BU):
        self.BU = BU
        self.comp = {}
        self.compute_compositions()

    def compute_compositions(self):
        # the fuel is divided into 3 equal segments.
        divisions  = ['fuel1', 'fuel2', 'fuel3']
        with open(comp_loc, 'rb') as f:
            stored_data =  pickle.load(f, encoding='latin1')
        # burnup steps
        BU = stored_data['tot']['BURNUP']
        # total volume of the fuel element
        total_volume = stored_data['tot']['VOLUME'][0]  # total volume of fuel element
        # total atomic density at  each time step
        tot_Adens = stored_data['tot']['ADENS']['total']
        # interploating to get the total atomic density at a given burnup value
        total_Adens_interp = interpolate.interp1d(BU, tot_Adens)(self.BU)*total_volume
        Adens ={}
        noise = ['total', 'lostdata', '];']
        for seg in divisions:
            # atomic densities at each segment
            Adens[seg] = stored_data[seg]['ADENS']
        Adens_accum = 0
        for  isotope in Adens['fuel1'].keys():
            if isotope not in noise:
                value = 0
                for seg in divisions:
                        # volume of the current segment
                        v = stored_data[seg]['VOLUME'][0]
                        # interpolating to get the atomic density of the istopee at this segment
                        # at given burnup
                        fun = interpolate.interp1d(BU, Adens[seg][isotope])
                        # getting number of atoms in this segment
                        data = fun(self.BU)*v
                        # adding the number of atoms of an isotope in each segment
                        value += data
                        Adens_accum += data
                # compute the atomic density of the isotope over the total fuel element
                self.comp[isotope] =  value/total_volume
                # computing the atomic denstity of H, ZR assuming H-Zr ratio 1.7
                #self.comp['H1'] = (total_Adens_interp - Adens_accum)*1.7/2.7/total_volume
                #self.comp['Zr90'] = (total_Adens_interp - Adens_accum)*1/2.7/total_volume
                # the next line was used when the recent data was chosen
                self.comp['H1'] = (total_Adens_interp - Adens_accum)/total_volume
