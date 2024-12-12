import numpy as np
import os
# remove the following when attenuation function completed
import random
import xraylib as xrl

THICK_MACRO = 'THICK'
MAT_MACRO = 'MAT'
TESTING = False

#Using local density definitions until package/library found with longer list
#in g/cm^3 
DENSITY = {'Si': 2.33, 'TiO2': 4.23,  'InSb': 5.78} 

def binary_representation_rev(num, bits):
    # need to reverse order so 1st bit in list is for first filter
    return format(num, '0'+str(bits)+'b')[::-1]

def get_densities(materials):
    
    densities = dict.fromkeys(materials)
    for material in list(densities):
        try:
            matdb = xrl.CompoundParser(material)
        except ValueError as err:
            print(f"{material} not found in xraylib.")
        else:
            if matdb['nAtomsAll'] == 1:
                density = xrl.ElementDensity(matdb['Elements'][0])
            else:
                if material in list(DENSITY):
                    density = DENSITY[material]
                else:
                    raise ValueError(f"{material} not found in DENSITY keys.")
            densities[material]=density
        
    return densities

def find_nearest_above(my_array, target):
    diff = my_array - target
    mask = np.ma.less(diff, 0)
    # We need to mask the negative differences and zero
    # since we are looking for values above
    if np.all(mask):
        return -1 # returns None if target is greater than any value
    masked_diff = np.ma.masked_array(diff, mask)
    return masked_diff.argmin()

class filterBlock():
    
    def __init__(self, filter_id = 'FL-1'):
        '''
        
        '''
        #create necessary inputs/output; initialize to zero/empty
        print(80*'#')
        print('Initializing pyDevice filter block object')
        self.materials = None
        self.thicknesses = None
        self.densities = None
        self.lookupTable = []
        self.culledTable = []
        self.sorted_index = []
        self.filter_attenuations = []
       
        self.filter_id = filter_id
        print(f'Filter ID: {self.filter_id}')
        
        self.outMask = 0
        self.inMask = 0
        self.config = 0
        self.attenuation = 1
        self.attenuation_actual = 1
        self.culledSize = 0
        self.configIndex = 0
        self.culledIndex = 0
        self.configIndexSorted = 0
        self.culledIndexSorted = 0
        self.energy = 0     #in keV
        
        self.verbose = False 
        
        self.num_filters = 0
        
        print(80*'#')   


    def setupLookupTable(self, subs_file, n_filters, energy = 8.0):
        '''
        lookup table created after IOC startup (after filter materials and 
        thicknesses are set
        '''
        print(80*'#')
        print('Setting up filter control...')
        
        self.num_filters = n_filters
        
        self.energy = energy
        
        #read in substitutions file
        try:
            subsFile = open(subs_file,"r")
        except:
            raise RuntimeError(f"Substiution file ({subsFile}) not found.")
        subsFileContent = subsFile.readlines()
        subsFile.close()
        
        macros = subsFileContent[2].replace('{','').replace('}','').replace(',','').split()
        filter_properties = {key: [] for key in macros} # dictionary of lists
        for i in range(self.num_filters):
            try:
                xx = subsFileContent[3+i].replace('{','').replace('}','').replace(',','').replace('"','').split()
                filter_properties[macros[0]].append(xx[0])
                filter_properties[macros[1]].append(xx[1])
                filter_properties[macros[2]].append(xx[2])
                filter_properties[macros[3]].append(xx[3])
            except:
                raise RuntimeError(f"Number of filters ({self.num_filters}) doesn't match substitution file")
        
        self.materials = []
        self.thicknesses = []

            
        # get materials from filter properties dictionary-list
        print('Getting filter materials...')
        if MAT_MACRO in macros:
            self.materials = filter_properties[MAT_MACRO]
            print('Filter material read in.\n')
        else:
            raise RuntimeError(f"Material macro ({MAT_MACRO}) not found in substituion file")
        
        # get densities from local definition (for compounds) or from xraylib (for elements)
        densities = get_densities(self.materials)
        self.densities = [densities[material] for material in self.materials]

        # get thicknesses from filter properties dictionary-list
        print('Getting filter thicknesses...')
        if THICK_MACRO in macros:
            self.thicknesses = [float(i) for i in filter_properties[THICK_MACRO]]
            print('Filter thicknesses read in.\n')
        else:
            raise RuntimeError(f"Thickness macro ({THICK_MACRO}) not found in substituion file")

        print('Calculating lookup table...')
        self.calc_lookup_table()
        print('Lookup table calculation complete.\n')
        
        print('Filter control setup complete.')
        print(80*'#')


    def calc_lookup_table(self):
        '''
        
        '''
        if TESTING:
            # For testing, using a shuffled version of the reciprocal of array whose values equal their index position + 1
            print('Using test array')
            self.lookupTable = np.reciprocal((np.arange(2**self.num_filters)+1).astype(float))
            np.random.shuffle(self.lookupTable)
        else:
            print('Using calculated attenuations')
            self.filter_attenuations = np.empty(self.num_filters)
            
            num_configs = 2**self.num_filters
            att_log = np.empty(num_configs)
            self.lookupTable = np.empty(num_configs)
            
            for i in range(self.num_filters):
                mass_thickness = self.densities[i]*(self.thicknesses[i]*1e-4)
                print(f'Getting attenuation for {self.materials[i]} at energy of {self.energy} of keV')
                self.filter_attenuations[i] = xrl.CS_Total_CP(self.materials[i], self.energy)*mass_thickness
                print(f'Filter: {i+1}, Mass_thickness: {mass_thickness}, log(attenuation): {self.filter_attenuations[i]}')
            for j in range(num_configs):
                binary_representation = binary_representation_rev(j, self.num_filters) 
                choice = [int(bit) for bit in binary_representation]
                att_log[j] = sum(self.filter_attenuations*choice)
            
            self.lookupTable = np.exp(att_log)
        
        self.cull_lookup_table()
            
    def cull_lookup_table(self):
        '''
        Culls the lookup table based on filters that are locked and or disabled
        '''
        self.culledSize = 2**(self.num_filters - (self.outMask | self.inMask).bit_count())
        if self.verbose: print(f'Operating spaced now at {self.culledSize} configurations')
        if self.verbose: print(f'Culling table with in mask {self.inMask} and out mask {self.outMask}')
        
        self.culledConfigs = np.empty(self.culledSize, dtype=int)
        self.culledTable = np.empty(self.culledSize)
        j = 0
        for i in range(2**self.num_filters):
            if ((i & self.outMask == 0) and (i & self.inMask == self.inMask)):
                self.culledConfigs[j]=i
                self.culledTable[j] = self.lookupTable[i]
                j += 1
                
        self.sort_lookup_table()
        
    def sort_lookup_table(self):
        '''
        
        '''
        if self.verbose: print(f'Sorting culled lookup table of length {len(self.culledTable)}')        
        self.sorted_index = np.argsort(self.culledTable)
        
    def calc_config_atten(self, filter_config, energy):
        '''
        Given filter set and incoming energy, calculate the attenuation
        '''
        att_log_total = 0
        
        binary_representation = binary_representation_rev(filter_config, self.num_filters)
        choice = [int(bit) for bit in binary_representation]
        for i in range(self.num_filters):
            if choice[i]:
                mass_thickness = self.densities[i]*(self.thicknesses[i]*1e-4)
                att_log_total += xrl.CS_Total_CP(self.materials[i], energy)*mass_thickness
        
        return np.exp(att_log_total)
        
    def updateEnergy(self, energy):
        '''
        Energy variable sent from IOC as a string
        '''
         
        self.energy = float(energy)
        
        if self.materials is None or self.thicknesses is None:
            raise RuntimeError("Substitution file hasn't been loaded. Can't calc lookup table")
        else:       
            self.calc_lookup_table()
            pydev.iointr(self.filter_id+'_new_energy', self.energy)
        self.setAttenActual()
        self.updateAttenRBVs()
        
    def updateAtten(self, attenuation):
        '''
        Attenuation variable sent from IOC as a string
        '''
        self.attenuation = float(attenuation)
        if self.verbose: print(f'Updating attenuation to {self.attenuation}')
        self.findFilterConfig()

    def findFilterConfig(self):
        ''' 
        User selects attenuation, this function finds nearest acheivable attenuation
        '''
        # Code to search lookup table for nearest attenuation to desired
        if self.verbose: print(f'Searching for config closest to {self.attenuation}')
        nearest_above_idx = find_nearest_above(self.culledTable, self.attenuation)
        if self.verbose: print(f'Nearest above index: {nearest_above_idx}')
        if nearest_above_idx == -1:
            # Desired attenuation is larger than possible at current energy and filter set
            self.culledIndex = np.argmax(self.culledTable)  
            if self.verbose: print(f'Desired attenuation greater than maximum possible; setting to {self.culledIndex}')
        else:
            self.culledIndex = nearest_above_idx

        #self.culledIndex = np.argmin(np.abs(self.culledTable - self.attenuation))
        if self.verbose: print(f'Config index found at {self.culledIndex}')

        self.culledIndexSorted = self.sorted_index.tolist().index(self.culledIndex)
        if self.verbose: print(f'Sorted config index found at {self.culledIndexSorted}')

        # Update PVs
        self.setAttenActual()
        self.updateFilterConfigPV()
        self.updateFilterRBV()
        self.updateAttenRBVs()
              
    def getPreviewAtten(self, sortedIndex):
        '''
        
        '''
        attenuation_preview = self.culledTable[self.sorted_index[int(sortedIndex)]]
        if self.verbose: print(f'Preview attenuation for {sortedIndex} is {attenuation_preview}')
        pydev.iointr(self.filter_id+'_new_preview', attenuation_preview)
    
    def updateIndex(self, sortedIndex):
        '''
        User has updated desired sorted index
        '''
        self.culledIndexSorted = int(sortedIndex)
        self.culledIndex = self.sorted_index[self.culledIndexSorted]
        self.setAttenActual()
        self.updateFilterConfigPV()
        self.updateFilterRBV()
        self.updateAttenRBVs()
        
    def updateConfig(self, config_BW):
        '''
        When user manually changes filters, this gets attenuation and displays it
        along with updated RBVs but it doesn't set the config PV
        '''
        self.config = int(config_BW)
        if self.verbose: print(f'User has set a filter, new base-2 config: {self.config}')
        
        self.culledIndex = (np.where(self.culledConfigs == self.config))[0][0]
        if self.verbose: print(f'User has set a filter, new culled index: {self.culledIndex}')
        self.culledIndexSorted = self.sorted_index.tolist().index(self.culledIndex)
        if self.verbose: print(f'User has set a filter, new sorted culled index: {self.culledIndexSorted}')
      
        self.setAttenActual()
        self.updateFilterRBV()
        self.updateAttenRBVs()

    def setInMask(self, inMask):
        '''
        update mask for filters that are locked in
        '''
        self.inMask = int(inMask)
        self.cull_lookup_table()
        if self.verbose: print(f'Converting culled index via in Mask')
        self.convertCulledIndex()
        if self.verbose: print(f'Updating filter RBV via in Mask')
        self.updateFilterRBV()
        if self.verbose: print(f'Setting in mask RBV to {self.inMask}')
        pydev.iointr(self.filter_id+'_new_inMask', int(self.inMask))

    def setOutMask(self, outMask):
        '''
        update mask for filters that must remain out (either disabled or locked)
        '''
        self.outMask = int(outMask)
        self.cull_lookup_table()
        if self.verbose: print(f'Converting culled index via out Mask')
        self.convertCulledIndex()
        if self.verbose: print(f'Updating filter RBV via out Mask')
        self.updateFilterRBV()
        if self.verbose: print(f'Setting out mask RBV to {self.outMask}')
        pydev.iointr(self.filter_id+'_new_outMask', int(self.outMask))

    def convertCulledIndex(self):
        '''
        When available configs change, need to update index so that tweaks 
        continue to work
        '''
        if self.verbose: print('Converting ...')
        self.culledIndex = (np.where(self.culledConfigs == self.config))[0][0]
        if self.verbose: print(f'Culled index is {self.culledIndex}')
        self.culledIndexSorted = self.sorted_index.tolist().index(self.culledIndex)
        if self.verbose: print(f'Sorted culled index is {self.culledIndexSorted}')
        
    def setAttenActual(self):
        '''
        Find the RBVs for the actual attenuation and the values of the next 
        up/down as a preview for the tweak buttons
        '''
        self.attenuation_actual = self.culledTable[self.culledIndex] 
        next_up = find_nearest_above(self.culledTable, self.culledTable[self.culledIndex])
        culledTable_sortedIndices = np.argsort(self.culledTable)
        sorted_position = np.where(culledTable_sortedIndices == self.culledIndex)[0][0]
        if self.verbose: print(f'Sorted position is {sorted_position}')
        self.next_up_atten = self.culledTable[culledTable_sortedIndices[min(sorted_position + 1,len(culledTable_sortedIndices)-1)]]
        self.next_dn_atten = self.culledTable[culledTable_sortedIndices[max(0,sorted_position - 1)]]
        if self.verbose: print(f'Next up: {self.next_up_atten} and next down: {self.next_dn_atten}')
        
    def updateFilterConfigPV(self):
        '''
        Update the filter config PV with the configuration associated 
        with the desired attenuation
        '''
        self.config = self.culledConfigs[self.culledIndex]
        pydev.iointr(self.filter_id+'_new_filters', int(self.config))
 
    def updateFilterRBV(self):
        '''
        Update the culled (masked) index RBV
        '''
        pydev.iointr(self.filter_id+'_new_index', int(self.culledIndexSorted))

    def updateAttenRBVs(self):
        '''
        Updated RBVs for various attenuation-related PVs
        '''
        pydev.iointr(self.filter_id+'_new_atten', self.attenuation_actual)
        pydev.iointr(self.filter_id+'_atten_up', self.next_up_atten)
        pydev.iointr(self.filter_id+'_atten_dn', self.next_dn_atten)
        self.attenuation_2E_actual = self.calc_config_atten(self.culledConfigs[self.culledIndex], self.energy*2)
        pydev.iointr(self.filter_id+'_new_atten_2E', self.attenuation_2E_actual)
        self.attenuation_3E_actual = self.calc_config_atten(self.culledConfigs[self.culledIndex], self.energy*3)
        pydev.iointr(self.filter_id+'_new_atten_3E', self.attenuation_3E_actual)
        
    def updateVerbosity(self, verbosity):
        '''
        Turn on minor printing
        '''
        print(f'Verbosity set to {verbosity}')
        self.verbose = verbosity
