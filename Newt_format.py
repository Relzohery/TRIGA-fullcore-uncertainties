import numpy as np
import json


Isotopes_of_interest = ['U234','U235','U236','U238','Np236','Np237','Pu238','Pu239','Pu240',
'Pu241','Pu242','Pu243','Pu244','Am241','Am242m','Am243','Cm242','Cm243','Cm244','Cs133','Cs134',
'Cs135','Cs137','Gd154','Gd155','Gd156','Tc99','Mo95','Kr83','Zr90','Zr91','Zr95','Zr93','Zr96',
'Zr94','Nb93','Mo97','Mo99','Ru101','Pd105','Ag109','I127','I129','Xe131','Nd143','Nd145','Sm147',
'Sm149','Sm150','Sm151','Sm152','Eu153','Eu154','Eu155','Ru103','Xe135','Pm147','Rh103','Ce144',
'Cd113','Eu151']

class NewtInput( ):
    '''' This class builds a full triga 2_D model used in the Triton depletion sequence (T-DEPL) '''
    def __init__(self, config, Isotopes_of_interest, core_comp = 'core_comp.json', power = 24.2653, burn = 2250, nlibs = 100, lib ='v7-56', Inputname = 'Input.inp'):
        self.config = config
        self.power = power
        self.core_comp = core_comp
        self.burn = burn
        self.nlibs = nlibs
        self.lib = lib
        self.celldata = True
        self.mats = {}
        self.locs = {}
        self.Input = Inputname
        self.make_materials()
        self.get_elements()
        self.IsotopesOfInterest= Isotopes_of_interest
        self.dimensions = {'zr_rod': 0.2286  , 'fuel_radius':  1.8161
             ,'clad_radius':1.8415,
             'control_radius':  1.5875,
             'control_clad':1.4224,
             'core_radius': 21.9075,
             'Al_outradius':22.86,
             'reflect_outradius': 30.0101 ,
             'Boxside':  32 }
        self.ringRadii = {'A':(0, 1) ,
                          'B': (4.05384, 1) ,
                          'C': (7.9806799999999996, 2) ,
                          'D': (11.945620000000002, 3),
                          'E': (15.91564, 4),
                          'F': (19.8882, 5)}

    def make_header(self):
        s = '=T-DEPL\n' # sequence name
        s += 'Triga full core 2-D\n' #title
        s += self.lib # library used
        s += '\n'
        return s


    def get_elements(self):
        # read configration file and skip the header
        lines = open(self.config,'r').readlines()[1: ]
        for line in lines:
            loc, ID = line.split()
            # make dictionary contains fuel Ids and their locations
            self.locs[loc] = ID

    def findLocation(self, location):
        ''' converts locations in core to x,y coordinates '''
        angles = np.linspace(0, 2 * np.pi, self.ringRadii[location[:1]][1] * 6, endpoint=False)
        radius = self.ringRadii[location[:1]][0]
        angle = angles[int(location[1:]) - 1]
        x = radius * np.cos(angle)
        y = radius * np.sin(angle)
        return (x, y)

    def loctoindex(self, loc):
        ''' converts  a location in core to index to be used as uint number and fuel materia numbers'''

        ringIndex = {l: i + 1 for i, l in enumerate('ABCDEF')}
        return '{}{:>02}'.format(ringIndex[loc[0]], loc[1:])

    def get_material_ids(self):
        ''' collects ids of fuel materais which are the same as units numbers '''
        mat = []
        for loc, Id in self.locs.items():
            if Id.isdigit():
                mat.append(self.loctoindex(loc))
        return mat

    def make_UnitsandHoles(self):
        ''' This makes units containing different elements, and also returns their representations as holes to be used in the global unit'''
        units = ''
        holes = ''
        controlRods = ['pulse', 'reg', 'shim', 'safety']

        for loc, ID in self.locs.items():
            ix = self.loctoindex(loc)
            coord = self.findLocation(loc)
            holes += 'hole {} origin x= {} y= {}\n'.format(ix, *coord)
            if ID.isdigit():
                # fule element region
                clad = str(self.mats['SS'][0][0])
                Zr = str(self.mats['Zr'][0][0])
                if self.celldata:
                    # if cell data is used, then all material inside the fuel unit should have unique ids
                    # so the material dict has to updated with those unique ids
                    self.mats['water'][0].append(str(self.mats['water'][0][0]) + str(ix))
                    clad = clad + str(ix)
                    self.mats['SS'][0].append(clad)
                    Zr = Zr + str(ix)
                    self.mats['Zr'][0].append(Zr)
                units += 'unit {}\n cylinder 1  {}\n'.format(ix, self.dimensions['clad_radius'])
                units+= 'cylinder 2 {}\ncylinder 3 {}\n'.format(self.dimensions['fuel_radius'], self.dimensions['zr_rod'])
                units+='media {}  1 3\nmedia {} 1 2 -3\nmedia {} 1 1 -2 -3\n'.format(Zr, ix,clad)
            elif ID.lower()  in controlRods:
                units += 'unit {}\n cylinder 1  {}\n'.format(ix, self.dimensions['control_clad'])
                units += 'cylinder 2  {}\n'.format(self.dimensions['control_radius'])
                units += 'media {} 1 2 \n'.format(self.mats['Boron'][0][0])
                units += 'media {}  1 -2 1\n'.format(self.mats['SS'][0][0])
            elif ID.lower() == 'graphite':
                units += 'unit {}\n cylinder 1  {}\n'.format(ix,self. dimensions['clad_radius'])
                units += '\n cylinder 2  {}\n'.format(self.dimensions['fuel_radius'])
                units += 'media {} 1 2 \n'.format(self.mats['graphite'][0][0])
                units += 'media {} 1 -2 1\n'.format(self.mats['SS'][0][0])
            else:
                units += 'unit {}\n cylinder 1  {}\n'.format(ix, self.dimensions['clad_radius'])
                units += 'media 1 1 1\n'
            units+= 'boundary 1 2 2\n'
        return units, holes


    def make_globalUnit(self):
        ''' This returns the global unit that surronds every thing in the reactor '''
        s = 'global unit 1\n'
        s += 'cylinder  1 {}\n'.format(self.dimensions['core_radius'])
        s += 'cylinder 2  {}\n'.format (self.dimensions['Al_outradius'])
        s += 'cylinder 3 {}\n'.format(self.dimensions['reflect_outradius'])
        s += 'cuboid 4 4p{}\n'.format(self.dimensions['Boxside'])
        s += self.make_UnitsandHoles()[1]
        s += 'media {}  1 1\n'.format(self.mats['water'][0][0])
        s += 'media {}  1 2 -1\n'.format(self.mats['Al'][0][0])
        s += 'media {} 1 3 -2 -1\n'.format(self.mats['graphite'][0][0])
        s += 'media {} 1  4 -3 -2 -1\n'.format( self.mats['water'][0][0])
        s += 'boundary 4  100 100\n'
        return s

    def read_geometry(self):
        # geometry block
        s  = ''
        s1 = self.make_UnitsandHoles()[0]
        s2 = self.make_globalUnit()
        s+= 'read geom\n'
        s +='{} {}'.format(s1, s2)
        s +='end geom\n'
        return s


    def read_model(self):

       ''' This returns the newt module  part that is placed in the middle of the depletion sequence'''

       s = 'read model\n2-D triga fuel pin\n'
       s += self.read_materials()
       s += self.read_geometry()
       s += self.read_bounds()
       s += 'end model\nend'

       return s

    def read_materials(self):
        # materiasl block
        s = 'read materials\n'
        # all materials other than fuel
        for mat_Ids in self.mats.values():
            for Id in mat_Ids[0]:
                s+= 'mix={}  end\n'.format(Id)
        # fuel material
        for loc, ID in self.locs.items():
            if ID.isdigit():
                Ix = self.loctoindex(loc)
                s +='mix={}  end\n'.format(Ix)
        s += 'end materials\n'
        return s

    def read_celldata(self):
        ''' This returns cell data block for cross section self sheliding, annular square pich cell is used '''

        s = 'read celldata\n'
        for Id in self.get_material_ids():
            clad = str(self.mats['SS'][0][0]) + str(Id)
            water = str(self.mats['water'][0][0]) + str(Id)
            Zr = str(self.mats['Zr'][0][0]) + str(Id)
            s+='Latticecell Asquarepitch pitch=4.1466 {} Fuelr 1.816 {} cladr=1.8415 {} IMODr=0.2286 {} end\n'.format(water,Id, clad, Zr)
        s += 'end celldata\n'
        return s

    def make_comp(self):
        ''' This returns the compositions block which is part of XSproc module '''

        s = 'read comp \n'
        for mat, mat_data in self.mats.items():
            for Id in mat_data[0]:
                for iso, comp in mat_data[1].items():
                    if mat == 'water':
                        s += iso + ' ' + str(Id) + '  DEN =' + str(comp) + ' end\n'
                    else:
                        s += iso + ' ' + str(Id) + ' ' + '0 '+ str(comp) + ' end\n'
        # fuel compositons  for each fuel element in this configuration are stored in json file,
        with open(self.core_comp) as core_data:
            core_comp = json.load(core_data)
        # fuel composiotns
        for loc, ID in self.locs.items():
            if ID.isdigit():
                Ix = self.loctoindex(loc) #fuel_material_Id
                element_comp = core_comp[ID]
                for isotope, Adens in element_comp.items():
                    # thermal scattering effect of hydrogen
                    if isotope.lower() == 'h1':
                        iso = 'h-zrh2'
                    # thermal scattering effect of hydrogen
                    elif isotope.lower() == 'zr90':
                        iso = 'zr90-zr5h8'
                    elif isotope.lower() =='am242m':
                        iso = 'Am-242m'
                    else:
                        iso = self.separate_str(isotope)
                    s += str(iso) + ' ' + str(Ix) + ' 0 '+ str(Adens)+ '  end\n'
        s += 'end comp\n'
        return s

    def read_depl(self):
        ''' This assigns material to be depleted '''

        s = ' read depletion \n '
        for loc, ID in self.locs.items():
            if ID.isdigit():
                Ix = self.loctoindex(loc)
                s += Ix + ' '
        s += '\nend depletion\n'
        return s

    def read_burndata(self):
        s = 'read burndata\n'
        s += 'power={}  burn={}  nlib={}  end\n'.format(self.power, self.burn, self.nlibs)
        s += 'end burndata\n'
        return s

    def read_opus(self):
        ''' This returns a block specifying the isotopes and the units desired '''
        s = 'read opus\n'
        s += '  title="Isotopic compositions in barn-cm - days"\n '
        s += 'units= Atom\n'
        s += 'symnuc= '
        s +=' '.join(self.IsotopesOfInterest)
        s+= ' end\n'
        s+='time=days\n'
        s+= 'end opus\n'
        return s

    def read_bounds(self):
         ''' This return a block that defines the boundary conditions '''
         s = 'read bounds\n'
         # using vacume boundary conditions in x and y directions
         s += 'yfc=vac   xfc=vac\n'
         s += 'end bounds\n'
         return s



    def make_input(self):
         ''' returns  a string with the whole T-depl depletion sequence '''
         s =''
         s+= self.make_header()
         # fill the matrial dict with fuel materials, cell matreials.
         self.make_UnitsandHoles()
         s+= self.make_comp()
         if self.celldata:
             s+=self.read_celldata()
         s+= self.read_depl()
         s+=self.read_burndata()
         s+=self.read_opus()
         s+=self.read_model()
         with open (self.Input, 'w') as f:
             f.write(s)



    def separate_str(self, string):
         digit =''
         s = ''
         for i in string:
             if i.isdigit():
                 digit += i
             else:
                s+= i
         return s+'-'+digit

    def make_materials(self):
        self.mats['water'] = ([1],{'water': 1.0})
        self.mats['Zr'] = ([2],{'Zr-90':  0.022077113304662624 ,
                                'Zr-91':  0.004814484184223803,
                                'Zr-92': 0.007359037768220872,
                                'Zr-94':  0.007457730403013338,
                                'Zr-96': 0.0012014755539952447 })

        self.mats['SS'] = ([3],{'C-12' :  0.00015932375780334336,
                'C-13' : 1.7232024749780397e-06 ,
                'Si-28': 0.0007939428979468555,
                'Si-29': 4.033291561628897e-05,
                'Si-30': 2.6618863412073752e-05,
                'P-31': 3.59083985829626e-05,
                'S-32': 2.1488343611263962e-05,
                'S-33': 1.6966267721284315e-07,
                'S-34':  9.614218375394447e-07,
                'S-36':  2.262169029504576e-09,
                'Cr-50': 0.000767775950402404,
                'Cr-52': 0.01480579496162648 ,
                'Cr-53': 0.0016788582979915406,
                'Mn-55': 0.0008802154068516285,
                'Fe-53': 0.0035516593020504528,
                'Fe-56': 0.055753455534702685,
                'Fe-57': 0.0012875904296056305,
                'Fe-58': 0.00017135464896120227,
                'Ni-58': 0.00518818463488483,
                'Ni-60': 0.00199846887613416 ,
                'Ni-61': 8.687238957805472e-05,
                'Ni-62': 0.00027699481284358066,
                'Ni-64': 7.0532850736459e-05})

        self.mats['Boron'] = ([5], {'B-11' : 0.10574621316977216 ,
                                    'B-10': 0.026271531112090706})


        self.mats['graphite'] = ([4], {'C-12':   0.11061493696317623,
                                        'C-13': 0.001196381103311418,
                                        'B-10': 2.4719626320659187e-08,
                                        'B-11': 9.949960142134675e-08 })
        self.mats['Al'] = ([6], {'Al-27': 0.060262620055259605})

if __name__== '__main__':
    Input = NewtInput('Config.txt', Isotopes_of_interest)
    Input.make_input()
    #with open('./Triga_2D.inp', 'w') as f:
        #f.write(s)
