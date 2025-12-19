import numpy as np
import tomllib
import xraylib
from collections import Counter
from transfocator_calcs import materials_to_deltas, materials_to_linear_attenuation
from transfocator_calcs import find_levels, get_densities
from transfocator_calcs import calc_1x_lu_table, calc_2x_lu_table, calc_kb_lu_table
from transfocator_calcs import calc_2xCRL_focus
from transfocator_calcs import SYSTEM_TYPE

OE_MACRO = 'OE'
MAT_MACRO = 'MAT'
NLENS_MACRO = 'NUMLENS'
RADIUS_MACRO = 'RADIUS'
LOC_MACRO = 'LOC'
THICKERR_MACRO = 'THICKERR'
DIM_MACRO = 'DIM'

modes = ['flat','round']

SYSTEM_TYPE_NAMES = {'1x': '1', '2x': '2', 'KB': '3'}

'''
Config variables

Beam Properties
        energy      : energy in keV
        L_und       : undulator length in m
        sigmaH_e    : Sigma electron source size in H direction in m
        sigmaV_e    : Sigma electron source size in V direction in m
        sigmaHp     : Sigma electron divergence in H direction in rad
        sigmaVp_e   : Sigma electron divergence in V direction in rad
Beamline properties
        d_StoL1 : Source-to-CRL1 distance, in m
        d_StoL2 : Source-to-CRL2 distance, in m
        d_Stof  : Source-to-focus distance, in m
CRL properties
        d_min   : Minimum thickness at the apex in m
        stacks  : number stacks in systems
KB properties
        ...     : ...
        ...     : ...
'''


DEFAULT_CONFIG = {'beam':{'energy': 15, 'L_und': 4.7, 
                          'round': {'sigmaH_e': 12.296e-6, 'sigmaV_e': 8.263e-6, 'sigmaHp_e': 2.336e-6, 'sigmaVp_e': 3.474e-6},
                          'flat':  {'sigmaH_e': 14.466e-6, 'sigmaV_e': 3.075e-6, 'sigmaHp_e': 2.749e-6, 'sigmaVp_e': 1.293e-6}},
                  'beamline': {'d_StoL': [51.9, 62.1], 'd_Stof': [55.0, 66.2]},
                  'crl':{'label': ['B', 'C'], 'stacks': [10, 10], 'd_min': [3.0e-5, 3.0e-5]},
                  'kb':{'label' : 'KB', 'KBH_L': 180.0e-3, 'KBH_q': 380.0e-3, 'KB_theta': 2.5e-3,
                        'KBV_L': 300.0e-3, 'KBV_q': 640.0e-3, 'KBH_p_limit': 1.0, 
                        'KBV_p_limit': 1.0 },
                  'sample': {'label': ['C', 'D']},
                  'init':{'sysType': '1x', 'initConfig': ['B']}}

DEFAULT_1X_CONFIG = {'sysType': SYSTEM_TYPE.singleCRL,
                     'CRLs': DEFAULT_CONFIG['crl']['label'][0],
                     'KBs': None,
                     'Sample': DEFAULT_CONFIG['sample']['label'][0]}

DEFAULT_2X_CONFIG = {'sysType': SYSTEM_TYPE.doubleCRL,
                     'CRLs': [DEFAULT_CONFIG['crl']['label'][0],DEFAULT_CONFIG['crl']['label'][1]],
                     'KBs': None,
                     'Sample': DEFAULT_CONFIG['sample']['label'][1]}
                     
DEFAULT_1XKB_CONFIG = {'sysType': SYSTEM_TYPE.CRLandKB,
                     'CRLs': DEFAULT_CONFIG['crl']['label'][0],
                     'KBs': DEFAULT_CONFIG['kb']['label'],
                     'Sample': DEFAULT_CONFIG['sample']['label'][1]}

def find_key(input_dict, target_value):
    '''
    Find key by value in a dictionary   
    '''
    key = [k for k, v in input_dict.items() if v == target_value]
    if key:
        return key[0]
   
def separate_by_oe(property_list, oe_list, desired_oe):
    '''
    Description:
        Lens properties are read in from substitutions file but not separated by 
        which transfocator they belong to.  This functions separates by optical 
        elemente (i.e. transfocator)
    
    Parameters:
        property_list   : list containing a property of all lenses
        oe_list         : list containing transfocator assignement of all lenses
        desired_oe      : which oe does user wnat properties for
        
    Returns:
        list of desired_oe's property values
    
    '''
    return [prop for (prop, oe) in zip(property_list, oe_list) if oe == desired_oe]
        


class focusingSystem():
    def __init__(self, crl_setup = None, beam_config = DEFAULT_CONFIG['beam'],
                 beamline_config = DEFAULT_CONFIG['beamline'],
                 crl_configs = DEFAULT_CONFIG['crl'], 
                 kb_config = DEFAULT_CONFIG['kb'],
                 sam_config = DEFAULT_CONFIG['sample'], 
                 init_config = DEFAULT_CONFIG['init']):

        '''
        Description:
            Focusing system object -- either single CRL, double CRL or single
            CRL + KB.
        
        Parameters:
           crl_setup        : ... 
                              Default: None 
           beam_config      : ... 
                              Default: DEFAULT_CONFIG['beam']
           beamline_config  : ... 
                              Default: DEFAULT_CONFIG['beamline']
           crl_configs      : ... 
                              Default: DEFAULT_CONFIG['crl'] 
           kb_config        : ... 
                              Default: DEFAULT_CONFIG['kb']
           init_config      : initial configuration of elements to be used
                              Default: DEFAULT_CONFIG['init']
            
       
        '''

        self.verbose = True

        self.elements = []
        self.toml_crls = []
        
        if crl_setup is None:
            beam = beam_config
            beamline = beamline_config
            crl = crl_configs
            kb = kb_config
            sample = sam_config
            init = init_config
            self.toml_file = 'DEFAULT'
        else:
            self.toml_file = crl_setup
            with open(crl_setup, "rb") as f:
                config = tomllib.load(f)
            beam = config['beam']
            beamline = config['beamline']
            #toml uses dictionary of lists; need to convert to list of dictionaries
            crls = config['crl']
            crl = []
            for i, label in enumerate(crls['labels']):
                self.toml_crls.append(label)
                crl.append({'label': label, 'stacks':crls['stacks'][i], 'd_min': crls['d_min']})

            # Add kbs if listed in toml
            if "kb" in config:
                kb = config['kb']
            else:
                kb = None
            
            # Add sample stations if listed in toml
            if "sample" in config:
                sample = config['sample']
                 
            init = config['init']


        # Create default configs
        self.single_config = DEFAULT_1X_CONFIG
        if len(crl) > 1:
            self.double_config = DEFAULT_2X_CONFIG
        if kb is not None:
            self.crlkb_config = DEFAULT_1XKB_CONFIG
        
        # Setup beam properties
        self.mode = modes[0]
        self.beam = {}
        self.setupSource(beam)
        
        # Setup element properties -- converting list of elements to dictionary  
        # keyed by labels
        self.crl = {}
        self.setupCRL(crl)

        # KB setup
        # KB systems need extra output info object distances for each mirror
        if kb is not None:
            self.kb = {}
            self.KB_ol = {'KBH_p_list': [],
                            'KBV_p_list': []}        
            self.setupKB(kb)        

        # Add possible sample stations
        self.sampleSTNs = sample['labels']

        # Setup beamline position of elements (CRLs and Samples) 
        # needs to come after self.crl and self.sampleSTNs sorted out
        self.bl = {}
        self.bl['d_StoL'] = {}
        self.bl['L_offset'] = {}
        self.bl['d_Stof'] = {}
        self.bl['f_offset'] = {}
        self.setupBeamline(beamline, crls['labels'])

            
        # Initialize slit sizes to 1 m (wide open)
        self.slits = {}
        self.setupSlits(kb = kb)

        # System type determined by init entry in config dictionary
        match init['sysType']:
            case SYSTEM_TYPE.singleCRL:
                self.curr_config = self.single_config
            case SYSTEM_TYPE.doubleCRL:
                self.curr_config = self.double_config
            case SYSTEM_TYPE.CRLandKB:
                self.curr_config = self.crlkb_config
        self.curr_config['CRLs'] = init['initConfig']
        
        self.updateElements()
         #<----------------------------------------------------------------------          
        # Are these needed at initialization        
        self.focalSize = 0 # get value from an ao (desired focal length)
        self.lenses = 0 # sets integer (2^10) whose binary representation indicates which lenses are in or out
        
        #---------------------------------------------------------------------->
        #initialize dictionary for crl indices of current state
        self.indexSorted = {'1':0, '2':0}
        if self.curr_config['sysType'] is SYSTEM_TYPE.doubleCRL:
            self.index = {'1':0,'2':0}
        else:
            self.index = {'1':0}
                        
        self.lookupTable = []
        self.q_list = []
        self.dq_list = []
        
        self.thickerr_flag = True
    
    def updateElements(self):
        self.elements = [elem for elem in self.curr_config['CRLs'] + [self.curr_config['KBs']] if elem is not None]
                
    def setupSource(self, beam_properties):
        '''
        Beam properties can have entries for the following
        
        energy      : energy in keV
        L_und       : undulator length in m
        sigmaH_e    : Sigma electron source size in H direction in m
        sigmaV_e    : Sigma electron source size in V direction in m
        sigmaHp_e   : Sigma electron divergence in H direction in rad
        sigmaVp_e   : Sigma electron divergence in V direction in rad
        
        
        mode        : two values: flat or round
        '''
        
        self.setEnergy(beam_properties['energy'])
        self.L_und = beam_properties['L_und']
        self.sigmaH_e = beam_properties['sigmaH_e']
        self.sigmaV_e = beam_properties['sigmaV_e']
        self.sigmaHp_e = beam_properties['sigmaHp_e']
        self.sigmaVp_e = beam_properties['sigmaVp_e']
                
        self.setupSourceEnergyDependent()

    def setEnergy(self, energy):
        '''
        Sets various forms of energy
        '''
        if energy > 0.0001:
            self.energy = float(energy)
            self.energy_eV = self.energy*1000.0  # Energy in keV
            self.wl = 1239.84 / (self.energy_eV * 10**9)    #Wavelength in m
            if self.verbose: print(f'Setting energy to {self.energy} keV')

    def setupSourceEnergyDependent(self, mode='flat'):
        '''
        Sets various energy dependent source parameters. Called whenever energy 
        is updated
        '''    
        _mode = (mode == modes[1])
        self.beam['sigmaH'] =  (self.sigmaH_e[modes[_mode]]**2 +  self.wl*self.L_und/2/np.pi/np.pi)**0.5
        self.beam['sigmaV'] =  (self.sigmaV_e[modes[_mode]]**2 +  self.wl*self.L_und/2/np.pi/np.pi)**0.5
        self.beam['sigmaHp'] = (self.sigmaHp_e[modes[_mode]]**2 + self.wl/self.L_und/2)**0.5
        self.beam['sigmaVp'] = (self.sigmaVp_e[modes[_mode]]**2 + self.wl/self.L_und/2)**0.5

    def setupBeamline(self, beamline_properties, crl_labels):
        '''
        Beamline properties can contain entries for the following
        
        d_StoL  : Source-to-CRLs distance, in m
        d_Stof  : Source-to-samples distance, in m
        '''            
        crl_distances = beamline_properties['d_StoL'] if isinstance(beamline_properties['d_StoL'], list) else [beamline_properties['d_StoL']]
        for crl_label, crl_dist in zip(crl_labels, crl_distances):
            self.bl['d_StoL'][crl_label] = crl_dist
            self.bl['L_offset'][crl_label] = 0
 
        sam_distances =  beamline_properties['d_Stof'] if isinstance(beamline_properties['d_Stof'], list) else [beamline_properties['d_Stof']]
        for sam_label, sam_dist in zip(self.sampleSTNs, sam_distances):
            self.bl['d_Stof'][sam_label] = sam_dist
            self.bl['f_offset'][sam_label] = 0

#       if self.sysType is singleCRLandKB # KB doesn't have location???
#           self.bl['d_StoKB'] = beamline_properties['d_StoKB']
            
    def setupCRL(self, crl):
        '''
        Looks through crl (list of transforcators) for entries for the following
        
        d_min   : Minimum thickness at the apex in m
        stacks  : number of stacks in system
        '''
        
        for tf in crl:
            self.crl[tf['label']]= {'d_min': tf['d_min'], 'stacks': tf['stacks']}
        
    
    def setupKB(self, kb):
        '''
        Looks through kb for kb properties

        KBH_L       : KBH length
        KBH_q       : KBH q
        KB_theta    : KB mirror angle
        KBV_L       : KBV length
        KBV_q       : KBV q
        KBH_p_limit : Minimum p that KBH can achieve, so abs(KBH_p) > KBH_p_min
        KBV_p_limit : Minimum p that KBV can achieve, so abs(KBV_p) > KBV_p_min
        '''
        
        self.kb['KBH_L']        = kb['KBH_L']
        self.kb['KBH_q']        = kb['KBH_q']
        self.kb['KB_theta']     = kb['KB_theta']
        self.kb['KBV_L']        = kb['KBV_L']
        self.kb['KBV_q']        = kb['KBV_q']
        self.kb['KBH_p_limit']  = kb['KBH_p_limit']
        self.kb['KBV_p_limit']  = kb['KBV_p_limit']
        
    def setupSlits(self, kb = None):
        '''
        Initializes slit sizes to 1 m (very open) 
        Assumes slits for all optical elements       
        '''
        for elem in self.toml_crls:
            self.slits[elem] = {'hor':1,'vert':1}
        if kb is not None:
            self.slits['KB'] = {'hor':1,'vert':1}
        
        
    def updateSlitSize(self, size, oe, slit):
        '''
        Slit size updates are propagated to CRL object from EPICS.  The beam
        size lookup table is then recalculated.
        '''
        
        self.slits[oe][slit] = float(size)
        if self.verbose: print(f"{oe} {slit} slit is set to {self.slits[oe][slit]}")
                             
    def updateSlitSizeRBV(self, oe, slit):
        '''
        Update proper slit size
        '''
        
        oes = oe if isinstance(oe, list) else [oe]

        for element in oes: 
            intr_string = 'updated_slitSize_'+element+'_'+slit
            pydev.iointr(intr_string, float(self.slits[element][slit]))
            if self.verbose: print(f"{oe} {slit} slit size RBV udpated to {self.slits[element][slit]}")
        
    def parseSubsFile(self, subs_file):
        '''
        Description:
        
        Parameters:
            ...             : ...
            
        Returns:
            ...             : ...
        
        '''
        #read in substitutions file
        try:
            subsFile = open(subs_file,"r")
        except:
            raise RuntimeError(f"Substitution file ({subsFile}) not found.")
        # Remove empty lines and comments
        subsFileContent = [line for line in subsFile.readlines() if (line.strip() and not line.startswith('#'))]
        subsFile.close()
        
        macros = subsFileContent[2].replace('{','').replace('}','').replace(',','').split()
        lens_properties = {key: [] for key in macros} # dictionary of lists

        i = 3
        while not subsFileContent[i].startswith('}'):
            try:
                xx = subsFileContent[i].replace('{','').replace('}','').replace(',','').replace('"','').split()
                lens_properties[macros[0]].append(xx[0])
                lens_properties[macros[1]].append(xx[1])
                lens_properties[macros[2]].append(xx[2])
                lens_properties[macros[3]].append(xx[3])
                lens_properties[macros[4]].append(xx[4])
                lens_properties[macros[5]].append(xx[5])
                lens_properties[macros[6]].append(xx[6])
                lens_properties[macros[7]].append(xx[7])
                lens_properties[macros[8]].append(xx[8])
            except:
                raise RuntimeError(f"""
                Substitutions file ({subs_file}) format error. 
                Last line attempted: {i}
                Line content: {xx}
                """)
            i += 1
        
        self.numlens = []
        self.radius = []
        self.materials = []
        self.lens_loc = []
        self.lens_thickerr = []

            
        # get number of lens for each lens stack from lens properties dictionary-list
        print('Getting OE assignments...')
        if OE_MACRO in macros:
#            self.oe_num = np.array([int(i) for i in lens_properties[OE_MACRO]])
            self.oe = lens_properties[OE_MACRO]
            print('OE assignments read in.\n')
        else:
            raise RuntimeError(f"OE assignemnt macro ({OE_MACRO}) not found in substituion file")
        
        oe_cnt = Counter(self.oe)
        self.subs_crls = list(oe_cnt.keys())
        # list of number of stacks for each crl
        self.num_stacks = [oe_cnt[crl] for crl in self.subs_crls]
            
        # get number of lens for each lens stack from lens properties dictionary-list
        print('Getting lens counts...')
        if NLENS_MACRO in macros:
            self.numlens = np.array([int(i) for i in lens_properties[NLENS_MACRO]])
            print('Number of lens read in.\n')
        else:
            raise RuntimeError(f"Number of lenses macro ({NLENS_MACRO}) not found in substituion file")

        # get radii for each lens from lens properties dictionary-list
        print('Getting lens\' radii...')
        if RADIUS_MACRO in macros:
            self.radius = np.array([float(i) for i in lens_properties[RADIUS_MACRO]])
            print('Radius of lenses read in.\n')
        else:
            raise RuntimeError(f"Radius macro ({RADIUS_MACRO}) not found in substituion file")

        # get materials from lens properties dictionary-list
        print('Getting lens materials...')
        if MAT_MACRO in macros:
            self.materials = lens_properties[MAT_MACRO]
            print('Lens material read in.\n')
        else:
            raise RuntimeError(f"Material macro ({MAT_MACRO}) not found in substituion file")
        
        # get densities from local definition (for compounds) or from xraylib (for elements)
        densities = get_densities(self.materials)
        self.densities = np.array([densities[material] for material in self.materials])

        # get location of each lens from lens properties dictionary-list
        print('Getting lens\' locations...')
        if LOC_MACRO in macros:
            self.lens_locations = np.array([float(l) for l in lens_properties[LOC_MACRO]])
#            self.lens_locations = np.array([float(l)*self.crl[str(self.oe_num[i])]['stack_d'] for i,l in enumerate(lens_properties[LOC_MACRO])])
            print('Location of lenses read in.\n')
        else:
            raise RuntimeError(f"Location macro ({LOC_MACRO}) not found in substituion file")


        # get dimensionality of each lens from lens properties dictionary-list
        print('Getting lens\' dimensionality...')
        if DIM_MACRO in macros:
            self.lens_dim = lens_properties[DIM_MACRO]
            print('Dimensionality of lenses read in.\n')
        else:
            raise RuntimeError(f"Dimensionality macro ({DIM_MACRO}) not found in substituion file")


        # get thicknesses errprfrom lens properties dictionary-list
        print('Getting lens thickness error...')
        if THICKERR_MACRO in macros:
            self.lens_thickerr = np.array([float(i) for i in lens_properties[THICKERR_MACRO]])
            print('Lens thickness errors read in.\n')
        else:
            raise RuntimeError(f"Thickness errors macro ({THICKERR_MACRO}) not found in substituion file")
            
        # Consistency check 1 
        # elements: subs file vs toml file (toml file read in during init)
        if self.toml_crls == self.subs_crls:
            self.crl_labels = self.toml_crls
        else:
            raise ValueError(f"""
            Substitution and toml files don't agree on optical elements.
            toml file ({self.toml_file}) has elements: {self.toml_crls}
            subs file ({subs_file}) has elements: {self.subs_crls}                      
            """)
        # Consistency check 2
        # number of stacks per element
        toml_stacks = [self.crl[s]['stacks'] for s in self.crl_labels]
        if self.num_stacks != toml_stacks:
            raise ValueError(f"""
            Substitution and toml files don't agree on number of stacks in each CRL:
            toml file ({self.toml_file}) has {toml_stacks} set of stacks
            subs file ({subs_file}) has {self.num_stacks} set of stacks          
            """)

    def setupLookupTable(self, subs_file):
        '''
        Description:
        
            -lookup table created after IOC startup
            -called directly by ioc startup file
            -Note: energy and slit size are updated before this is called but table 
        calculation is disabled during their updates
        
        Parameters
            subs_file   : EPICS substitutions file
            n_stacks    : list of number of stacks in each transfocator system
            
        '''
        print(80*'#')
        print('Setting up lens control...')
        
        self.parseSubsFile(subs_file)

        # Dictionary of config chosen for each element
        self.config = {}

        # Dictionary of total possible configs for each element
        self.configs = {}
        for c, n in zip(self.crl_labels, self.num_stacks):
            self.configs[c] = np.arange(2**n)
 #       for i, n in enumerate(self.num_stacks): self.configs[self.crl[i]['label']] = np.arange(2**n)
                
        self.lens_count = {}
        self.radii = {}
        self.mat = {}
        self.lens_loc = {}
        self.thickerr = {}
        for crl in self.crl_labels:
            self.lens_count[crl] = separate_by_oe(self.numlens, self.oe, crl)
            self.radii[crl] = separate_by_oe(self.radius, self.oe, crl)
            self.mat[crl] = separate_by_oe(self.materials, self.oe, crl)
            self.lens_loc[crl]  = separate_by_oe(self.lens_locations, self.oe, crl)
            self.thickerr[crl]  = separate_by_oe(self.lens_thickerr, self.oe, crl)


        print(f" Lens count: {self.lens_count}")
                 
        print('Constructing lookup table...')
        self.construct_lookup_table()
        print('Lookup table calculation complete.\n')
        
        print('Transfocator control setup complete.')
        print(80*'#')
        
    def construct_lookup_table(self):
        '''
        Description:
            Constructs lookup table for three types of systems:
            -single transfocator
            -double transfocator
            -transfocator + KB mirror
            Updates IOC waveforms and rBS with output of table constructions
            Should be called after beam energy or slits size changes
        '''
        if self.verbose: print(f"Constructing lookup table for {self.curr_config['sysType']} system")
        match self.curr_config['sysType']:
            case SYSTEM_TYPE.singleCRL:

                crl = self.curr_config['CRLs'][0]
                sam = self.curr_config['Sample']
                
                bl_subset = {'d_StoL1': bl['d_StoL'][crl],
                             'L1_offset': bl['L_offset'][crl],
                             'd_Stof': bl['d_Stof'][sam],
                             'f_offset': bl['f_offset'][sam]}
                                
                results_dict = calc_1x_lu_table(self.configs[crl].size(), self.radii[crl], self.mat[crl], 
                                           self.energy, self.wl, self.lens_count[crl], 
                                           self.lens_loc[crl], self.beam, bl_subset, 
                                           self.crl[crl], self.slits[crl]['hor'], 
                                           self.slits[crl]['vert'], self.thickerr[crl], 
                                           flag_HE = self.thickerr_flag, verbose = self.verbose)

            case SYSTEM_TYPE.doubleCRL:
                crl1 = self.curr_config['CRLs'][0]
                crl2 = self.curr_config['CRLs'][1]
                sam = self.curr_config['Sample']
 

                bl_subset = {'d_StoL1': bl['d_StoL'][crl1],
                             'L1_offset': bl['L_offset'][crl1],
                             'd_StoL2': bl['d_StoL'][crl2],
                             'L2_offset': bl['L_offset'][crl2],
                             'd_Stof': bl['d_Stof'][sam],
                             'f_offset': bl['f_offset'][sam]}
                    

                results_dict = calc_2x_lu_table(self.configs[crl1].size(), self.radii[crl1], 
                                            self.mat[crl1], self.radii[crl2], 
                                            self.mat[crl2], self.energy, self.wl,
                                            self.lens_count, self.lens_loc[crl1],
                                            self.lens_loc[crl2], self.beam, bl_subset,
                                            self.crl, self.slits, self.thickerr[crl1], 
                                            self.thickerr[crl2], flag_HE = self.thickerr_flag,
                                            verbose = self.verbose)
                
                
                self.index1to2_sorted = results_dict['invf2_indices']           
            case SYSTEM_TYPE.CRLandKB:
                crl = self.curr_config['CRLs'][0]
                sam = self.curr_config['Sample']

                bl_subset = {'d_StoL1': bl['d_StoL'][crl],
                             'L1_offset': bl['L_offset'][crl],
                             'd_Stof': bl['d_Stof'][sam],
                             'f_offset': bl['f_offset'][sam]}

                results_dict = calc_kb_lu_table(self.configs[crl].size(), self.radii[crl],
                                            self.mat[crl], self.energy, self.wl,
                                            self.lens_count[crl], self.lens_loc[crl], 
                                            self.beam, bl_subset, self.crl, self.kb, 
                                            self.slits, self.thickerr[crl], 
                                            flag_HE = self.thickerr_flag,
                                            verbose = self.verbose)


                self.KB_ol = {'KBH_p_list': results_dict['KBH_p_list'], 
                            'KBV_p_list': results_dict['KBV_p_list']}
                            
        self.lookupTable = results_dict['FWHM_atsample_list']
        self.sorted_invF_index = results_dict['invF_list_sort_indices']
        self.sorted_invF = results_dict['invF_list_sorted']                                                          
        self.q_list = results_dict['q_list']
        self.dq_list = results_dict['dq_list']
                                                                    
        self.updateEnergyRBV()
        self.updateModeRBV()

        self.updateSlitSizeRBV(self.elements, 'hor')
        self.updateSlitSizeRBV(self.elements, 'vert')

        self.updateLookupWaveform()
        self.updateInvFWaveform()
        self.updateLookupConfigs()  
        
        if self.curr_config['sysType'] == SYSTEM_TYPE.doubleCRL: 
            self.setFocalSizeActual(offTable = True)
        else:
            self.setFocalSizeActual(offTable = False)
        self.updateFocalSizeRBVs()

        self.updateQWaveforms()
        self.updateQdistances()
        
        if self.curr_config['sysType'] == SYSTEM_TYPE.CRLandKB:
            self.updateKBWaveforms()
            self.updateKBdistanceRBVs()

    def updateSysType(self, sysType):
        # Save current configuration
        match self.curr_config:
            case {'sysType': SYSTEM_TYPE.singleCRL, **rest}:
                self.single_config = self.curr_config
            case {'sysType': SYSTEM_TYPE.doubleCRL, **rest}:
                self.double_config = self.curr_config
            case {'sysType': SYSTEM_TYPE.CRLandKB, **rest}:
                self.crlkb_config = self.curr_config

        # Update current system type
        self.curr_config['sysType'] = sysType
    
        # Restore configuration for that system type
        match sysType:
            case SYSTEM_TYPE.singleCRL:
                self.curr_config = self.single_config
            case SYSTEM_TYPE.doubleCRL:
                self.curr_config = self.double_config
            case SYSTEM_TYPE.CRLandKB:
                self.curr_config = self.crlkb_config
        
        for i, crl in enumerate(self.curr_config['CRLs']):
            self.updateSystemRBV(i)       
        
        # TODO when KB to be fully integrated
#        if len(self.curr_config['KB']) > 0
#           kb_iointr_name = ...
#           pydev.iointr(kb_iointr_name, self.curr_config['KB'])

        self.updateElements()
        
        # update sample
        pydev.iointr('updated_sample', self.curr_config['Sample'])
        pydev.iointr('updated_sysType', self.curr_config['sysType'])
        
    def assignSystem(self, systemNum, oe):        
        self.curr_config['CRLs'][systemNum-1] = oe   

        self.updateElements()

        # Set readback
        self.updateSystemRBV(systemNum)       

    def updateSystemRBV(self, systemNum):          
        # Set readback       
        iointr_name = 'updated_system' + systemNum
        pydev.iointr(iointr_name, self.curr_config['CRLs'][systemNum-1])


    def assignSample(self, sampleSTN):
        
        if sampleSTN not in self.sampleSTNs:
            raise ValueError(f"""
            Requested sample station {sampleSTN} is not in toml file ({self.toml_file}) list.
            Pre-defined sample stations include: {self.sampleSTNS}. Check the sampleSTN DB file
            to confirm that it is consistent with the TOML file.        
            """)
        self.curr_config['Sample'] = sampleSTN   
  
        # Set readback       
        self.updateSampleRBV()

    def updateSampleRBV(self):
  
        # Set readback       
        iointr_name = 'updated_sample'
        pydev.iointr(iointr_name, self.curr_config['Sample'])

       
    def updateModeRBV(self):
        '''
        Description
            Updates beam mode readback PV.  To be called after lookup table calculated
        '''
        retval = (self.mode == modes[1])
        pydev.iointr('updated_mode', float(retval))

    def updateEnergyRBV(self):
        '''
        Description
            Updates energy readback PV.  To be called after lookup table calculated
        '''
        pydev.iointr('updated_E', float(self.energy))

    def updateInvFWaveform(self):
        '''
        Description
            Puts invF lists into waveform PVs after lookup table calculatione
        '''
        crl1_label = self.curr_config['CRLs'][0]
        crl2_label = self.curr_config['CRLs'][1]
        pydev.iointr('new_invFind_list_'+crl1_label, self.sorted_invF_index['1'].tolist())
        pydev.iointr('new_invF_list_'+crl1_label, self.sorted_invF['1'].tolist())
        if self.sorted_invF_index['2'] is not None:
            pydev.iointr('new_invFind_list_'+crl2_label, self.sorted_invF_index['2'].tolist())
            pydev.iointr('new_invF_list_'+crl2_label, self.sorted_invF['2'].tolist())
            
    def updateLookupConfigs(self):
        '''
        Description
            Puts lookup table config integers into waveform PV after lookup 
            table calculation
        '''        
        for label in self.curr_config['CRLs']:
            interrupt_str = 'new_configs_'+label
            pydev.iointr(interrupt_str, self.configs[label].tolist())
            

    def updateIndex(self, sortedIndex, oe):
        '''
        Description
            User has updated desired sorted index for either CRL1 or CRL2.  In 
            double CRL case, this update (either of CRL1 or CRL2**) moves the 
            system off the lookup table.

            If double CRL system's CRL1 index is moved, it will move on lookup 
            table, so CRL2 will be udpated as well

            **In practice, user should only move CRL2 in  the double case
            
        Parameters
            sortedIndex: string
                configuration index to set optical element to
            oe: string
                Label of optical element 
        '''
        if self.verbose: print(f'Setting {oe} to index {sortedIndex}')
        oe_num = str(self.curr_config['CRLs'].index(oe)+1)

        self.indexSorted[oe_num] = int(sortedIndex)
        if oe_num == '1':
            self.index['1'] = self.sorted_invF_index['1'][self.indexSorted['1']] 
            if self.curr_config['sysType'] == SYSTEM_TYPE.doubleCRL:
                self.indexSorted['2'] = self.index1to2_sorted[self.indexSorted['1']]
                self.index['2'] = self.sorted_invF_index['2'][self.indexSorted['2']]
        elif oe_num == '2':
            self.index['2'] = self.sorted_invF_index['2'][self.indexSorted['2']]

        # Update PVs
        if oe_num == '2': 
            self.setFocalSizeActual(offTable = True)
        else:
            self.setFocalSizeActual(offTable = False)

        self.updateLensConfigPV()
        self.updateLensRBV()
        self.updateFocalSizeRBVs()  
 
        self.updateQWaveforms()
        self.updateQdistances()
        
        # KB system need to "publish" p_h, p_v, and q1 
        if self.curr_config['sysType'] == SYSTEM_TYPE.CRLandKB:
            self.updateKBdistanceRBVs()

    def updateConfig(self, config_RBV, oe):
        '''
        Description
            Readback of lens configuration has changed: get focal size and 
            display it along with updated RBVs but don't set the config PV
            
        Parameters
            config_BW: string
                configuration as binary
            oe: string
                Label of optical element 
        '''

        if self.verbose: print(f'Getting focal size for {oe} set to {config_RBV}')
        oe_num = str(self.curr_config['CRLs'].index(oe)+1)

        self.index[oe_num] = config_RBV
        # Find the configuration in the 1/f sorted list
        self.indexSorted[oe_num] = self.sorted_invF_index[oe_num].tolist().index(self.index[oe_num])

        if self.curr_config['sysType'] == SYSTEM_TYPE.doubleCRL: 
            self.setFocalSizeActual(offTable = True)
        else:
            self.setFocalSizeActual(offTable = False)

        self.updateQdistances()
        self.updateLensRBV()
        self.updateFocalSizeRBVs()   

    def updateFsize(self, focalSize):
        '''
        Descriptoin:
            User updates desired focal size. Lookup table is traversed to find 
            nearest to desired.
            
        Parameters:
            focalSize: string
                Desired focal size (in m)
        '''
        # focalPoint variable sent from IOC as a string
        self.focalSize = float(focalSize)
        if self.verbose: print(f'Setting focal size to {self.focalSize}')
        self.find_config()

    def updateZpos(self, zOffset, elem, etype):
        '''
        Description
            User has updated the position of an optical element or the sample
            beamline parameters are updated
            
        Parameters
            zOffset: float
                distance from local zero for the element
            elem: string
                Label of element 
        '''
                
        if self.verbose: print(f'Shifting {elem} position by {zOffset}')
        if etype == 'crl':
            self.bl['L_offset'][elem]=float(zOffset)
            self.updateGeneralRBV(f'new_{elem}_Offset',zOffset)  
            self.updateGeneralRBV(f'new_{elem}_pos',float(zOffset)+self.bl['d_StoL'][elem])           
        elif etype == 'sam':
            self.bl['f_offset'][elem]=float(zOffset)
            self.updateGeneralRBV(f'new_sam{elem}PosOffset',zOffset)  
            self.updateGeneralRBV(f'new_sam{elem}Pos',float(zOffset)+self.bl['d_Stof'][elem])

                    
    def find_config(self):
        ''' 
        Description:
            User selected focal size, this function finds nearest acheivable 
            focal size from the lookup table
        '''
        # Code to search lookup table for nearest focal size to desired; note the
        # lookup table is already sorted by 1/f
        if self.verbose: print(f'Searching for config closest to {self.focalSize}')

        # simple approach
        # self.indexSorted = np.argmin(np.abs(self.lookupTable - self.focalSize))

        # XS approach -- can handle nan but in pydev application don't have a good
        # way to "transmit" errors (i.e. no solution found) to user.
        indices, _ = find_levels(self.lookupTable, self.focalSize, direction='forward')
        self.indexSorted['1'] = indices[0]
        if self.verbose: print(f"1/f-sorted config index found at {self.indexSorted['1']}")

        self.index['1'] = self.sorted_invF_index['1'][self.indexSorted['1']]
        if self.verbose: print(f"CRL 1 config index found at {self.index['1']}")
        
        if self.curr_config['sysType'] == SYSTEM_TYPE.doubleCRL:
            self.indexSorted['2'] = self.index1to2_sorted[self.indexSorted['1']]
            self.index['2'] = self.sorted_invF_index['2'][self.indexSorted['2']]
            if self.verbose: print(f"CRL 2 config index found at {self.index['2']}")
            

        # Update PVs
        if self.verbose: print(f'Updating RBVs')

        self.setFocalSizeActual(offTable = False)
        self.updateLensConfigPV()
        self.updateLensRBV()
        self.updateFocalSizeRBVs()     

        self.updateQWaveforms()
        self.updateQdistances()
        
        # KB system need to "publish" p_h, p_v, and q1 
        if self.curr_config['sysType'] == SYSTEM_TYPE.CRLandKB:
            self.updateKBdistanceRBVs()

    def setFocalSizeActual(self, offTable = False):
        '''
        Description:
            Gets the actual focal size for new system configuration
        Parameters:
            offTable: boolean
                For double CRL, when 2nd CRL is changed this parameter is True and 
                causes calculation of the focal size for the new CRL 2 + current 
                CRL 1 index
        '''

        if self.verbose: print(f'Setting actual focal size')
        if not offTable:
            self.focalSize_actual = self.lookupTable[self.indexSorted['1']] 
            self.q = self.q_list[self.indexSorted['1']]   
            self.dq = self.dq_list[self.indexSorted['1']]   
        else:
            crl1 = self.curr_config['CRLs'][0]
            crl2 = self.curr_config['CRLs'][1]
            sam = self.curr_config['Sample']

            bl_subset = {'d_StoL1': bl['d_StoL'][crl1],
                         'L1_offset': bl['L_offset'][crl1],
                         'd_StoL2': bl['d_StoL'][crl2],
                         'L2_offset': bl['L_offset'][crl2],
                         'd_Stof': bl['d_Stof'][sam],
                         'f_offset': bl['f_offset'][sam]}
        
            
            fsize, q2, dq2 = calc_2xCRL_focus(self.index[crl1], self.index[crl2], 
                                              self.radii[crl1], self.mat[crl1], 
                                              self.radii[crl2], self.mat[crl2], 
                                              self.energy, self.wl,
                                              self.lens_count, self.lens_loc[crl1],
                                              self.lens_loc[crl2],
                                              self.beam, bl_subset, self.crl, 
                                              self.slits,
                                              self.thickerr[crl1], 
                                              self.thickerr[crl2], 
                                              flag_HE = self.thickerr_flag,
                                              verbose = self.verbose)
            self.focalSize_actual = fsize
            self.q = q2
            self.dq = dq2
         
        
    def updateLensConfigPV(self):
        '''
        Description:
            Updates optical element config PVs for which stacks need to be in/out
        '''

        for i, crl_label in enumerate(self.curr_config['CRLs']):
            if self.verbose: print(f'Setting lens configuration PV for CRL {i+1}')
            self.config[crl_label] = self.index[str(i+1)]
            pydev.iointr('new_lenses_'+crl_label, int(self.config[crl_label]))            

    def updateLensRBV(self):
        '''
        Description;
            Updates optical elements config index PVs
        '''

        for i, crl_label in enumerate(self.curr_config['CRLs']):
            if self.verbose: print(f'Setting lens configuration index RBV for CRL {i+1}: {self.indexSorted[str(i+1)]} ')
            pydev.iointr('new_index_'+crl_label, int(self.indexSorted[str(i+1)]))            
                    
    def updateFocalSizeRBVs(self):
        '''
        Description:
            Updated focal size readback PV
        '''
        if self.verbose: print(f'Setting actual focal size to {self.focalSize_actual}')
        pydev.iointr('new_fSize', self.focalSize_actual)
 
    def updateQdistances(self):
    
        if self.verbose: 
            print(f'New image distance for last CRL:  {self.q}')
            print(f'New image distance to sample from last CRL:  {self.dq}')

        pydev.iointr('new_q', self.q)
        pydev.iointr('new_dq', self.dq)



    def updateKBdistanceRBVs(self):
        '''
        Description:
            Updated image/object distance readback PVs for KB system
        '''
        kbh_p = self.KB_ol['KBH_p_list'][self.indexSorted['1']]
        kbv_p = self.KB_ol['KBV_p_list'][self.indexSorted['1']]

        if self.verbose: 
            print(f'New object distance for horizontal KB:  {kbh_p}')
            print(f'New object distance for vertical KB:  {kbv_p}')

        pydev.iointr('KBH_p_list', kbh_p)
        pydev.iointr('KBV_p_list', kbv_p)

        
        
    def getPreviewFocalSize(self, sortedIndex):
        '''
        Description:
            Finds focal size for desired index
            
        Parameters:
            sortedIndex: string
                index user would like preview focal size
        '''
        fSize_preview = self.lookupTable[int(sortedIndex)]
        if self.verbose: print(f'Preview focal sizes for {sortedIndex} is {fSize_preview}')
        pydev.iointr('new_preview', fSize_preview)
    
    def setThickerrFlag(self, flag):
        '''
        Description:
            User has updated thickness error flag.
        
        Parameters:
            flag : string
                converted to boolean with 
                    True: each stack's thickness error is used in focal size calculation
                    False: stack thickness error NOT used in focal size calculation 
        '''
        self.thickerr_flag = int(flag)
        if self.verbose: print(f'Thickness Error Flag set to {flag}')
        self.updateThickerrFlagRBV()

    def updateThickerrFlagRBV(self):
        '''
        Description:
            Thickness error flag has been updated; readback PV is set
        '''
        if self.verbose: print(f'Thickness Error Flag RBV set to {self.thickerr_flag}')
        pydev.iointr('updated_thickerr_Flag', self.thickerr_flag)
 

    def updateMode(self, mode):
        '''
        Description
        
        '''
        
        if self.verbose: print(f'Updating beam mode from {self.mode} to {modes[mode]}')
        self.mode = modes[mode]
        self.setupSourceEnergyDependent(mode=self.mode)
        
    def updateE(self, energy):
        '''
        Description:
            Beam energy updates are propagated to CRL object from EPICS. The beam
            size lookup table is then recalculated.

        Parameters:
            energy: string
                beam energy in keV
        '''

        if energy > 0.0001:
            # Energy variable sent from IOC as a string
            self.setEnergy(energy)
            # Update beam properties that are dependent on energy
            self.setupSourceEnergyDependent()
        else:
            if self.verbose: print(f'Invalid energy setting: {energy} kev; staying at {self.energy} keV')
       
 
    def updateLookupWaveform(self):
        '''
        Description:
                Puts lookup table focal sizes into waveform PV
        '''
        pydev.iointr('new_lookupTable', self.lookupTable.tolist())

    def updateQWaveforms(self):
        '''
        Description:
                Puts object distances lists into waveform PVs
        '''
        pydev.iointr('new_q_list', self.q_list.tolist())
        pydev.iointr('new_dq_list', self.dq_list.tolist())


    def updateKBWaveforms(self):
        '''
        Description:
                Puts object/image distances lists for KB system into waveform PVs
        '''
        pydev.iointr('new_KBH_p_list', self.KB_ol['KBH_p_list'].tolist())
        pydev.iointr('new_KBV_p_list', self.KB_ol['KBV_p_list'].tolist())

    def updateGeneralRBV(self, interrupt_str, val):
        '''
        Description:
            Puts value to PV
            
        Parameters:
            interrupt_str   : string
                pydev interrupt label that PV is watching
            val             : scalar
                value to put into PV
        '''
        pydev.iointr(interrupt_str, val)


    def updateVerbosity(self, verbosity):
        '''
        Description:
            Turn on messages to iocConsole from python code
            
        Parameters:
            verbosity: string
            numerical (0 or 1, for now) value to later be used as boolean
        '''
        print(f'Verbosity set to {verbosity}')
        self.verbose = int(verbosity)

