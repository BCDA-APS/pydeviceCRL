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

SYSTEM_TYPE_NAMES = {'1x': '1', '2x', '2', 'KB': '3'}

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
                  'beamline': {'d_StoL1': 51.9, 'd_StoL2': 62.1, 'd_Stof': 66.2},
                  'crl':[{'stacks': 10, 'd_min': 3.0e-5}],
                  'kb':{'KBH_L': 180.0e-3, 'KBH_q': 380.0e-3, 'KB_theta': 2.5e-3,
                        'KBV_L': 300.0e-3, 'KBV_q': 640.0e-3, 'KBH_p_limit': 1.0, 
                        'KBV_p_limit': 1.0 }}

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
    # TODO: ms update
    def __init__(self, crl_setup = None, beam_config = DEFAULT_CONFIG['beam'],
                 beamline_config = DEFAULT_CONFIG['beamline'],
                 crl_configs = DEFAULT_CONFIG['crl'], 
                 kb_config = DEFAULT_CONFIG['kb'], sysType = SYSTEM_TYPE.singleCRL):

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
           sysType          : ... 
                              Default: SYSTEM_TYPE.singleCRL
            
       
        '''

        self.verbose = True

        # TODO: ms update
        self.elements = []
        self.n_elements = 0
        # TODO: ms update
        if crl_setup is None:
            beam = beam_config
            beamline = beamline_config
            crl = crl_configs
            kb = kb_config
            self.sysType = sysType
        else:
            # TODO: ms update
            with open(crl_setup, "rb") as f:
                config = tomllib.load(f)
            beam = config['beam']
            beamline = config['beamline']
            #check config for the crl, crl1, crl2, kb and use this to determine system type
            crl = []

            # TODO: ms update
            if "crl" in config:
                self.n_elements+=1
                crl.append(config['crl'])
                self.elements.append('1')
                if "kb" in config:
                    self.n_elements+=1
                    kb = config['kb']
                    self.elements.append('kb')
            # TODO: ms update
            if "crl1" in config:
                self.n_elements+=1 
                crl.append(config['crl1'])
                self.elements.append('1')
                if "crl2" in config:
                    self.n_elements+=1
                    crl.append(config['crl2'])
                    self.elements.append('2')
            if "init" in config:
               self.sysType = config['init']['sysType']
			   self.elements = config['init']['elems']
                
        
        # TODO: ms update
       # Setup beam properties
        self.mode = modes[0]
        self.beam = {}
        self.setupSource(beam)
        
        # TODO: ms update
        # Setup beamline position of elements
        self.bl = {}
        self.setupBeamline(beamline)

        # TODO: ms update
        # Setup element properties
        self.crl = {}
        self.setupCRL(crl)
        if self.sysType is SYSTEM_TYPE.CRLandKB:
            self.setupKB(kb)

        # TODO: ms update
        # Initialize slit sizes to 0
        self.slits = {}
        self.setupSlits()
      
        #<----------------------------------------------------------------------          
        # Are these needed at initialization        
        self.focalSize = 0 # get value from an ao (desired focal length)
        self.lenses = 0 # sets integer (2^10) whose binary representation indicates which lenses are in or out
        
        #---------------------------------------------------------------------->
        #initialize dictionary for crl indices of current state
        self.indexSorted = {'1':0, '2':0}
        if self.sysType is SYSTEM_TYPE.doubleCRL:
            self.index = {'1':0,'2':0}
        else:
            self.index = {'1':0}
        
        # KB systems need extra output info object distances for each mirror
        if self.sysType is SYSTEM_TYPE.CRLandKB:
            self.KB_ol = {'KBH_p_list': [],
                            'KBV_p_list': []}
                
        self.lookupTable = []
        self.q_list = []
        self.dq_list = []
        
        self.thickerr_flag = True

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

    # TODO: ms update
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

    # TODO: ms update
    def setupBeamline(self, beamline_properties, num=1):
        '''
        Beamline properties can contain entries for the following
        
        d_StoL1 : Source-to-CRL1 distance, in m
        d_StoL2 : Source-to-CRL2 distance, in m
        d_Stof  : Source-to-sample distance, in m
        '''            
        
        self.bl['d_StoL1'] = beamline_properties['d_StoL1']
        self.bl['L1_offset']=0
        self.bl['d_Stof'] = beamline_properties['d_Stof']
        self.bl['f_offset']=0
        if self.sysType is SYSTEM_TYPE.doubleCRL: 
            self.bl['d_StoL2'] = beamline_properties['d_StoL2']
            self.bl['L2_offset']=0
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
        
    
    # TODO: ms update
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
        
    # TODO: ms update
    def setupSlits(self):
        '''
        Initializes slit sizes to 1 m (very open)        
        '''
        self.slits['1'] = {'hor':1,'vert':1}            
        if self.sysType is SYSTEM_TYPE.doubleCRL:
            self.slits['2'] = {'hor':1,'vert':1}            
        if self.sysType is SYSTEM_TYPE.CRLandKB:
            self.slits['KB'] = {'hor':1,'vert':1}

        
    # TODO: ms update
    def updateSlitSize(self, size, oe, slit):
        '''
        Slit size updates are propagated to CRL object from EPICS.  The beam
        size lookup table is then recalculated.
        '''
        
        self.slits[oe][slit] = float(size)
        if self.verbose: print(f"{oe} {slit} slit is set to {self.slits[oe][slit]}")
                             
    # TODO: ms update
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
            raise RuntimeError(f"Substiution file ({subsFile}) not found.")
        # Remove empty lines and comments
        subsFileContent = [line for line in subsFile.readlines() if (line.strip() and not line.startswith('#'))]
        subsFile.close()
        
        macros = subsFileContent[2].replace('{','').replace('}','').replace(',','').split()
        lens_properties = {key: [] for key in macros} # dictionary of lists

        i = 3
        while subsFilesContent[i] != '}':
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
                raise RuntimeError(f"Substitutions file ({subsFile}) format error")
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
		
        oe_cnt = Counter(oe)
		self.elements = list(oe_cnt.keys())
        self.num_stacks = [oe_cnt[elem] for elem in self.elements]
            
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
            print('Location of lenses read in.\n')
        else:
            raise RuntimeError(f"Dimensionality macro ({DIM_MACRO}) not found in substituion file")


        # get thicknesses errprfrom lens properties dictionary-list
        print('Getting lens thickness error...')
        if THICKERR_MACRO in macros:
            self.lens_thickerr = np.array([float(i) for i in lens_properties[THICKERR_MACRO]])
            print('Lens thickness errors read in.\n')
        else:
            raise RuntimeError(f"Thickness errors macro ({THICKERR_MACRO}) not found in substituion file")

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
        
        # Since element configs are paired, the lookup table size will be equal
        # to minimum possible configs for one of the system's elements
        # Though shouldn't it be the # of US CRL configs?
        self.parseSubsFile(subs_file)
        self.num_configs = 2**(min(self.num_stacks))

        # Dictionary of config chosen for each element
        self.config = {}

        # Dictionary of total possible configs for each element
        self.configs = {}
        for i, n in enumerate(self.num_stacks): self.configs[self.elements[i]] = np.arange(2**n)
                
        self.lens_count = {}
        self.radii = {}
        self.mat = {}
        self.lens_loc = {}
        self.thickerr = {}
        for oe in self.oe:
			self.lens_count[oe] = separate_by_oe(self.numlens, self.oe, oe)
			self.radii[oe] = separate_by_oe(self.radius, self.oe, oe)
			self.mat[oe] = separate_by_oe(self.materials, self.oe, oe)
			self.lens_loc[oe]  = separate_by_oe(self.lens_locations, self.oe, oe)
			self.thickerr[oe]  = separate_by_oe(self.lens_thickerr, self.oe, oe)

        print(self.lens_count)
                 
        print('Constructing lookup table...')
        # TODO: ms update -- need to set which OEs to use where....maybe have
        # an initial setting?
        # By default setting to 1x, with 1st CRL as active optical element
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
    
        if self.sysType == SYSTEM_TYPE.singleCRL:
			elem1 = self.elements[0]
            results_dict = calc_1x_lu_table(self.configs[elem1], self.radii[elem1], self.mat[elem1], 
                                           self.energy, self.wl, self.lens_count[elem1], 
                                           self.lens_loc[elem1], self.beam, self.bl, 
                                           self.crl, self.slits[elem1]['hor'], 
                                           self.slits[elem1]['vert'], self.thickerr[elem1], 
                                           flag_HE = self.thickerr_flag, verbose = self.verbose)
                                                   
        elif self.sysType == SYSTEM_TYPE.doubleCRL:
			elem1 = self.elements[0]
			elem2 = self.elements[1]
 
            results_dict = calc_2x_lu_table(self.configs[elem1], self.radii[elem1], 
                                            self.mat[elem1], self.radii[elem2], 
                                            self.mat[elem2], self.energy, self.wl,
                                            self.lens_count, self.lens_loc[elem1],
                                            self.lens_loc[elem2], self.beam, self.bl,
                                            self.crl, self.slits, self.thickerr[elem1], 
                                            self.thickerr[elem2], flag_HE = self.thickerr_flag,
                                            verbose = self.verbose)
                                                   
            self.index1to2_sorted = results_dict['invf2_indices']
            
        elif self.sysType == SYSTEM_TYPE.CRLandKB:
			elem1 = self.elements[0]
            results_dict = calc_kb_lu_table(self.configs[elem1], self.radii[elem1],
                                            self.mat[elem1], self.energy, self.wl,
                                            self.lens_count[elem1], self.lens_loc[elem1], 
                                            self.beam, self.bl, self.crl, self.kb, 
                                            self.slits, self.thickerr[elem1], 
                                            flag_HE = self.thickerr_flag,
                                            verbose = self.verbose)

            self.KB_ol = {'KBH_p_list': results_dict['KBH_p_list'], 
                          'KBV_p_list': results_dict['KBV_p_list']}

            
        self.lookupTable = results_dict['FWHM_atsample_list']
        self.sorted_invF_index = results_dict['invF_list_sort_indices']
        self.sorted_invF = results_dict['invF_list_sorted']                                                          
        self.q_list = results_dict['q_list']
        self.dq_list = results_dict['dq_list']
                                                                    
        # TODO: ms update
        self.updateEnergyRBV()
        self.updateModeRBV()
        self.updateSlitSizeRBV(self.elements, 'hor')
        self.updateSlitSizeRBV(self.elements, 'vert')

        # TODO: ms update
        self.updateLookupWaveform()
        self.updateInvFWaveform()
        self.updateLookupConfigs()  
        
        # TODO: ms update
        if self.sysType == SYSTEM_TYPE.doubleCRL: 
            self.setFocalSizeActual(offTable = True)
        else:
            self.setFocalSizeActual(offTable = False)
        self.updateFocalSizeRBVs()

        # TODO: ms update
        self.updateQWaveforms()
        self.updateQdistances()
        
        # TODO: ms update
        if self.sysType == SYSTEM_TYPE.CRLandKB:
            self.updateKBWaveforms()
            self.updateKBdistanceRBVs()

    # TODO: ms update
    def updateSysType(self, sysType):
        self.sysType = SYSTEM_TYPE_NAMES[sysType]

        # set some defaults for element assignments?
        # or go to previous values?
        
        # Other object updates?
        # TODO
        # self.redoSetupLookupTable()
        # self.construct_lookup_table()
        

        pydev.iointr('updated_sysType', find_key(SYSTEM_TYPE_NAMES, self.sysType))
        
    def assignSystem(self, systemNum, oe):        
        self.elements[int(systemNum)-1] = oe
          
        # TODO: Other object updates?
        # self.redoSetupLookupTable()
        # self.construct_lookup_table()
  
        # Set readback       
        iointr_name = 'updated_system' + systemNum
        pydev.iointr(iointr_name, self.elements[int(systemNum)-1]))

       
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

        pydev.iointr('new_invFind_list_1', self.sorted_invF_index['1'].tolist())
        pydev.iointr('new_invF_list_1', self.sorted_invF['1'].tolist())
        if self.sorted_invF_index['2'] is not None:
            pydev.iointr('new_invFind_list_2', self.sorted_invF_index['2'].tolist())
            pydev.iointr('new_invF_list_2', self.sorted_invF['2'].tolist())
            
    def updateLookupConfigs(self):
        '''
        Description
            Puts lookup table config integers into waveform PV after lookup 
            table calculation
        '''
        
        pydev.iointr('new_configs_1', self.configs['1'].tolist())
        if self.sysType == SYSTEM_TYPE.doubleCRL:
            pydev.iointr('new_configs_2', self.configs['2'].tolist())
            

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
        self.indexSorted[oe] = int(sortedIndex)
        if oe == '1':
            self.index['1'] = self.sorted_invF_index['1'][self.indexSorted['1']] 
            if self.sysType == SYSTEM_TYPE.doubleCRL:
                self.indexSorted['2'] = self.index1to2_sorted[self.indexSorted['1']]
                self.index['2'] = self.sorted_invF_index['2'][self.indexSorted['2']]
        elif oe == '2':
            self.index['2'] = self.sorted_invF_index['2'][self.indexSorted['2']]
            
        # Update PVs
        if oe == '2': 
            self.setFocalSizeActual(offTable = True)
        else:
            self.setFocalSizeActual(offTable = False)

        self.updateLensConfigPV()
        self.updateLensRBV()
        self.updateFocalSizeRBVs()  
 
        self.updateQWaveforms()
        self.updateQdistances()
        
        # KB system need to "publish" p_h, p_v, and q1 
        if self.sysType == SYSTEM_TYPE.CRLandKB:
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

        self.index[oe] = config_RBV
        # Find the configuration in the 1/f sorted list
        self.indexSorted[oe] = self.sorted_invF_index[oe].tolist().index(self.index[oe])

        if self.sysType == SYSTEM_TYPE.doubleCRL: 
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

    def updateZpos(self, zOffset, elem):
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
        if elem == '1':
            self.bl['L1_offset']=float(zOffset)
            self.updateGeneralRBV('new_1_Offset',zOffset)  
            self.updateGeneralRBV('new_1_pos',float(zOffset)+self.bl['d_StoL1'])  
        elif elem == '2':
            self.bl['L2_offset']=float(zOffset)
            self.updateGeneralRBV('new_2_Offset',zOffset)  
            self.updateGeneralRBV('new_2_pos',float(zOffset)+self.bl['d_StoL2'])  
        elif elem == 'sam':
            self.bl['f_offset']=float(zOffset)
            self.updateGeneralRBV('new_samPosOffset',zOffset)  
            self.updateGeneralRBV('new_samPos',float(zOffset)+self.bl['d_Stof'])  
                    
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
        
        if self.sysType == SYSTEM_TYPE.doubleCRL:
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
        if self.sysType == SYSTEM_TYPE.CRLandKB:
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
            fsize, q2, dq2 = calc_2xCRL_focus(self.index['1'], self.index['2'], 
                                              self.radii['1'], self.mat['1'], 
                                              self.radii['2'], self.mat['2'], 
                                              self.energy, self.wl,
                                              self.lens_count, self.lens_loc['1'],
                                              self.lens_loc['2'],
                                              self.beam, self.bl, self.crl, 
                                              self.slits,
                                              self.thickerr['1'], 
                                              self.thickerr['2'], 
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
        if self.verbose: print(f'Setting lens configuration PV for CRL 1')
        self.config['1'] = self.index['1']
        pydev.iointr('new_lenses_1', int(self.config['1']))
        if self.sysType == SYSTEM_TYPE.doubleCRL:
            if self.verbose: print(f'Setting lens configuration PV for CRL 2')
            self.config['2'] = self.index['2']
            pydev.iointr('new_lenses_2', int(self.config['2']))
            

    def updateLensRBV(self):
        '''
        Description;
            Updates optical elements config index PVs
        '''
        if self.verbose: print(f"Setting lens configuration index RBV for CRL 1: {self.indexSorted['1']}")
        pydev.iointr('new_index_1', int(self.indexSorted['1']))
        if self.sysType == SYSTEM_TYPE.doubleCRL:
            if self.verbose: print(f"Setting lens configuration index RBV for CRL 2: {self.indexSorted['2']}")
            pydev.iointr('new_index_2', int(self.indexSorted['2']))
                    
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

