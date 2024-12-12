import numpy as np
import tomllib
import xraylib
from transfocator_calcs import lookup_diameter, materials_to_deltas, materials_to_linear_attenuation
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
        stack_d : Stack thickness in m
        stacks  : number stacks in systems
KB properties
        ...     : ...
        ...     : ...
'''
DEFAULT_CONFIG = {'beam':{'energy': 15, 'L_und': 4.7, 'sigmaH_e': 14.8e-6,
                          'sigmaV_e': 3.7e-6, 'sigmaHp_e': 2.8e-6, 'sigmaVp_e': 1.5e-6},
                  'beamline': {'d_StoL1': 51.9, 'd_StoL2': 62.1, 'd_Stof': 66.2},
                  'crl':[{'stacks': 10, 'stack_d': 50.0e-3, 'd_min': 3.0e-5}],
                  'kb':{'KBH_L': 180.0e-3, 'KBH_q': 380.0e-3, 'KB_theta': 2.5e-3,
                        'KBV_L': 300.0e-3, 'KBV_q': 640.0e-3, 'KBH_p_limit': 1.0, 
                        'KBV_p_limit': 1.0 }}
    
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

        self.elements = []
        self.n_elements = 0
        if crl_setup is None:
            beam = beam_config
            beamline = beamline_config
            crl = crl_configs
            kb = kb_config
            self.sysType = sysType
        else:
            with open(crl_setup, "rb") as f:
                config = tomllib.load(f)
            beam = config['beam']
            beamline = config['beamline']
            #check config for the crl, crl1, crl2, kb and use this to determine system type
            crl = []

            if "crl" in config:
                self.n_elements+=1
                crl.append(config['crl'])
                self.sysType = SYSTEM_TYPE.singleCRL
                self.elements.append('1')
                if "kb" in config:
                    self.n_elements+=1
                    kb = config['kb']
                    self.sysType = SYSTEM_TYPE.CRLandKB
                    self.elements.append('kb')
            if "crl1" in config:
                self.n_elements+=1 
                crl.append(config['crl1'])
                self.sysType = SYSTEM_TYPE.singleCRL
                self.elements.append('1')
                if "crl2" in config:
                    self.n_elements+=1
                    crl.append(config['crl2'])
                    self.sysType = SYSTEM_TYPE.doubleCRL
                    self.elements.append('2')
                
        
        # Setup beam properties
        self.beam = {}
        self.setupSource(beam)
        
        # Setup beamline position of elements
        self.bl = {}
        self.setupBeamline(beamline)

        # Setup element properties
        self.crl = {}
        self.setupCRL(crl)
        if self.sysType is SYSTEM_TYPE.CRLandKB:
            self.setupKB(kb)

        # Initialize slit sizes to 0
        self.slits = {}
        self.setupSlits()
      
        #<----------------------------------------------------------------------          
        # Are these needed at initialization
        
        #TODO -- any generalizations? Yes but how? Need to do by elements?
        self.energy = 0  # gets value from an ao (incoming beam energy)
        self.focalSize = 0 # get value from an ao (desired focal length)
        self.lenses = 0 # sets integer (2^10) whose binary representation indicates which lenses are in or out
        
        #---------------------------------------------------------------------->
        #initialize dictionary for crl indices of current state
        self.indexSorted = {'1':0, '2':0}
        if self.sysType is SYSTEM_TYPE.doubleCRL:
            self.index = {'1':0,'2':0}
        else:
            self.index = {'1':0}
        
                
        self.lookupTable = []
        
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

    def setupSourceEnergyDependent(self):
        '''
        Sets various energy dependent source parameters. Called whenever energy 
        is updated
        '''    
        self.beam['sigmaH'] =  (self.sigmaH_e**2 +  self.wl*self.L_und/2/np.pi/np.pi)**0.5
        self.beam['sigmaV'] =  (self.sigmaV_e**2 +  self.wl*self.L_und/2/np.pi/np.pi)**0.5
        self.beam['sigmaHp'] = (self.sigmaHp_e**2 + self.wl/self.L_und/2)**0.5
        self.beam['sigmaVp'] = (self.sigmaVp_e**2 + self.wl/self.L_und/2)**0.5

    def setupBeamline(self, beamline_properties, num=1):
        '''
        Beamline properties can contain entries for the following
        
        d_StoL1 : Source-to-CRL1 distance, in m
        d_StoL2 : Source-to-CRL2 distance, in m
        d_Stof  : Source-to-sample distance, in m
        '''            
        
        self.bl['d_StoL1'] = beamline_properties['d_StoL1']
        self.bl['d_Stof'] = beamline_properties['d_Stof']
        if self.sysType is SYSTEM_TYPE.doubleCRL: 
            self.bl['d_StoL2'] = beamline_properties['d_StoL2']
#       if self.sysType is singleCRLandKB # KB doesn't have location???
#           self.bl['d_StoKB'] = beamline_properties['d_StoKB']
            
    def setupCRL(self, crl):
        '''
        Looks through crl (list of transforcators) for entries for the following
        
        d_min   : Minimum thickness at the apex in m
        stack_d : Stack thickness in m
        stacks  : number of stacks in system
        '''
        
        for elem, tf in enumerate(crl):
            self.crl[str(elem+1)]= {'d_min': tf['d_min'], 'stack_d': tf['stack_d'], 'stacks': tf['stacks']}
        
    
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
        
    def setupSlits(self):
        '''
        Initializes slit sizes to 0        
        '''
        self.slits['1'] = {'hor':0,'vert':0}            
        if self.sysType is SYSTEM_TYPE.doubleCRL:
            self.slits['2'] = {'hor':0,'vert':0}            
        if self.sysType is SYSTEM_TYPE.CRLandKB:
            self.slits['KB'] = {'hor':0,'vert':0}

        
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
            raise RuntimeError(f"Substiution file ({subsFile}) not found.")
        # Remove empty lines and comments
        subsFileContent = [line for line in subsFile.readlines() if (line.strip() and not line.startswith('#'))]
        subsFile.close()
        
        macros = subsFileContent[2].replace('{','').replace('}','').replace(',','').split()
        lens_properties = {key: [] for key in macros} # dictionary of lists
        
        for i in range(self.total_stacks):
            try:
                xx = subsFileContent[3+i].replace('{','').replace('}','').replace(',','').replace('"','').split()
                lens_properties[macros[0]].append(xx[0])
                lens_properties[macros[1]].append(xx[1])
                lens_properties[macros[2]].append(xx[2])
                lens_properties[macros[3]].append(xx[3])
                lens_properties[macros[4]].append(xx[4])
                lens_properties[macros[5]].append(xx[5])
                lens_properties[macros[6]].append(xx[6])
                lens_properties[macros[7]].append(xx[7])
            except:
                raise RuntimeError(f"Number of lenses ({self.total_stacks}) doesn't match substitution file")
        
        self.numlens = []
        self.radius = []
        self.materials = []
        self.lens_loc = []
        self.lens_thickerr = []

            
        # get number of lens for each lens stack from lens properties dictionary-list
        print('Getting OE assignments...')
        if OE_MACRO in macros:
            self.oe_num = np.array([int(i) for i in lens_properties[OE_MACRO]])
            print('OE assignments read in.\n')
        else:
            raise RuntimeError(f"OE assignemnt macro ({OE_MACRO}) not found in substituion file")

            
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
            self.lens_locations = np.array([float(l)*self.crl[str(self.oe_num[i])]['stack_d'] for i,l in enumerate(lens_properties[LOC_MACRO])])
            print('Location of lenses read in.\n')
        else:
            raise RuntimeError(f"Location macro ({LOC_MACRO}) not found in substituion file")

        # get thicknesses errprfrom lens properties dictionary-list
        print('Getting lens thickness error...')
        if THICKERR_MACRO in macros:
            self.lens_thickerr = np.array([float(i) for i in lens_properties[THICKERR_MACRO]])
            print('Lens thickness errors read in.\n')
        else:
            raise RuntimeError(f"Thickness errors macro ({THICKERR_MACRO}) not found in substituion file")

    def setupLookupTable(self, subs_file, n_stacks):
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
        
        # convert number of lenses to list
        self.num_stacks = n_stacks if isinstance(n_stacks, list) else [n_stacks]

        self.total_stacks = sum(self.num_stacks)
        # Since element configs are paired, the lookup table size will be equal
        # to minimum possible configs for one of the system's elements
        self.num_configs = 2**(min(self.num_stacks))

        self.parseSubsFile(subs_file)

        # Dictionary of config chosen for each element
        self.config = {}

        # Dictionary of total possible configs for each element
        self.configs = {}
        for i, n in enumerate(self.num_stacks): self.configs[str(i+1)] = np.arange(2**n)
                
        self.lens_count = {}
        self.lens_count['1'] = separate_by_oe(self.numlens, self.oe_num, 1)
        self.lens_count['2'] = separate_by_oe(self.numlens, self.oe_num, 2)
        print(self.lens_count)
        
        self.radii = {}
        self.radii['1'] = separate_by_oe(self.radius, self.oe_num, 1)
        self.radii['2'] = separate_by_oe(self.radius, self.oe_num, 2)
        
        self.mat = {}
        self.mat['1'] = separate_by_oe(self.materials, self.oe_num, 1)
        self.mat['2'] = separate_by_oe(self.materials, self.oe_num, 2)
        
        self.lens_loc = {}
        self.lens_loc['1']  = separate_by_oe(self.lens_locations, self.oe_num, 1)
        self.lens_loc['2']  = separate_by_oe(self.lens_locations, self.oe_num, 2)
        
        self.thickerr = {}
        self.thickerr['1']  = separate_by_oe(self.lens_thickerr, self.oe_num, 1)
        self.thickerr['2']  = separate_by_oe(self.lens_thickerr, self.oe_num, 2)
        
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
    

        if self.sysType == SYSTEM_TYPE.singleCRL:
            arr_a, dict_b, dict_c = calc_1x_lu_table(self.num_configs, 
                                                   self.radii['1'], self.mat['1'], 
                                                   self.energy, self.wl,
                                                   self.lens_count['1'], self.lens_loc['1'], 
                                                   self.beam, self.bl, self.crl, 
                                                   self.slits['1']['hor'], self.slits['1']['vert'],
                                                   self.thickerr['1'], 
                                                   flag_HE = self.thickerr_flag,
                                                   verbose = self.verbose)
        elif self.sysType == SYSTEM_TYPE.doubleCRL: 
            arr_a, dict_b, dict_c, arr_d = calc_2x_lu_table(self.num_configs, 
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
            self.index1to2_sorted = arr_d
        elif self.sysType == SYSTEM_TYPE.CRLandKB:
            arr_a, dict_b, dict_c = calc_kb_lu_table(self.num_configs, 
                                                   self.radii['1'], self.mat['1'], 
                                                   self.energy, self.wl,
                                                   self.lens_count['1'], self.lens_loc['1'], 
                                                   self.beam, self.bl, self.crl,
                                                   self.kb, self.slits, 
                                                   self.thickerr['1'], 
                                                   flag_HE = self.thickerr_flag,
                                                   verbose = self.verbose)
            
        self.lookupTable = arr_a
        self.sorted_invF_index = dict_b
        self.sorted_invF = dict_c                                                            
                                                                    
        self.updateEnergyRBV()
        self.updateSlitSizeRBV(self.elements, 'hor')
        self.updateSlitSizeRBV(self.elements, 'vert')

        self.updateLookupWaveform()
        self.updateInvFWaveform()
        self.updateLookupConfigs()  
        
        if self.sysType == SYSTEM_TYPE.doubleCRL: 
            self.setFocalSizeActual(offTable = True)
        else:
            self.setFocalSizeActual(offTable = False)
        self.updateFocalSizeRBVs()       

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
        else:
            self.focalSize_actual = calc_2xCRL_focus(self.index['1'], self.index['2'], 
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
            if verbose: print(f'Invalid energy setting: {energy} kev; staying at {self.energy} keV')
       
 
    def updateLookupWaveform(self):
        '''
        Description:
                Puts lookup table focal sizes into waveform PV
        '''
        pydev.iointr('new_lookupTable', self.lookupTable.tolist())


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

