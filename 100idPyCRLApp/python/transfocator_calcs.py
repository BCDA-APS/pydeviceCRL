import xraylib
from scipy.optimize import root_scalar
import numpy as np
from enum import Enum

__all__ = """
    lookup_diameter 
    materials_to_deltas 
    materials_to_linear_attenuation 
    calc_lookup_table
    calc_single_lookup_table
    get_densities
""".split()

# Lookup table where each entry is a tuple (column1, column2)
Lens_diameter_table = [
    (50, 450.0),
    (100, 632.0),
    (200, 894.0),
    (300, 1095.0),
    (500, 1414.0),
    (1000, 2000.0),
    (1500, 2450.0),
]

# Convert the lookup table to a dictionary for faster lookup
Lens_diameter_dict = {int(col1): col2 for col1, col2 in Lens_diameter_table}

#Using local density definitions until package/library found with longer list
#in g/cm^3 
DENSITY = {'Si': 2.33, 'TiO2': 4.23,  'InSb': 5.78} 

class SYSTEM_TYPE(Enum):
    singleCRL = 1
    doubleCRL = 2 
    CRLandKB = 3
    
def get_densities(materials):
    '''
    Description:
        Gets densities for the lens materials used in transfocator
        
    Parameters:
        materials   : list of strings
            chemical symbols for lens materials
        
    Returns:
        densities: dictionary with keys set to the chemical symbol of lens materials
            and values set to density of lens materials
    
    '''
    
    densities = dict.fromkeys(materials)
    for material in list(densities):
        try:
            matdb = xraylib.CompoundParser(material)
        except ValueError as err:
            print(f"{material} not found in xraylib.")
        else:
            if matdb['nAtomsAll'] == 1:
                density = xraylib.ElementDensity(matdb['Elements'][0])
            else:
                if material in list(DENSITY):
                    density = DENSITY[material]
                else:
                    raise ValueError(f"{material} not found in DENSITY keys.")
            densities[material]=density
        
    return densities

def index_to_binary_list(index, length):
    """
    Converts an index number to its binary representation as a list of digits,
    and pads the list with zeros in front to achieve the desired length.
    
    Parameters:
        index (int): The index number to be converted.
        length (int): The desired length of the binary list.
    
    Returns:
        list: A list of digits representing the binary representation of the index.
    """
    # Convert the index to a binary string and remove the '0b' prefix
    binary_str = bin(index)[2:]
    
    # Pad the binary string with zeros in front to achieve the desired length
    #padded_binary_str = binary_str.zfill(length)
      
    # Reverse the binary string
    reversed_binary_str = binary_str[::-1]
    
    # Convert the reversed binary string to a list of integers
    binary_list = [int(digit) for digit in reversed_binary_str]
    
    # Pad the list with zeros at the end to achieve the desired length
    while len(binary_list) < length:
        binary_list.append(0)

    return binary_list

def binary_list_to_index(binary_list, length):
    """
    Converts a list of binary digits in reverse order to its integer representation,
    padding the list with zeros at the end to have a fixed number of elements.
    
    Parameters:
        binary_list (list): A list of digits representing the binary number in reverse order.
        length (int): The fixed number of elements the list should have.
    
    Returns:
        int: The integer representation of the binary number.
    """
    # Pad the list with zeros at the end to achieve the desired length
    while len(binary_list) < length:
        binary_list.append(0)
    
    # Convert the binary list to an integer
    index = 0
    for i, digit in enumerate(binary_list):
        index += digit * 2**i
        
    return index



def lookup_diameter(lens_radius):
    '''
    Description:
        Convert lens radius to a diameter (See XS for more)
    Parameters:
        lens_radius : float
            lens radius (in m) to be converted to diameter
    Returns:
        lens diameter   : float
            from lookup table, converted to m
    
    '''
    # Convert the input float to an integer
    input_int = int(round(lens_radius*1.0e6))
    return Lens_diameter_dict.get(input_int, (lens_radius*1.0e6)**0.5*63.222+ 0.73)/1.0e6

def materials_to_deltas(material_list, energy):
    """
    Description:
        Convert a list of material names to a list of delta values at a given energy.
    
    Parameters:
        material_list (list): A list of material names.
        energy (float): The energy in keV.
    
    Returns:
        list: A list of delta values for the given materials at the given energy.
    """
    # The list to store delta values
    delta_list = []

    # Iterate through each material in the input list
    for material in material_list:
        # Compute the delta value for the current material at the given energy
        Z = xraylib.SymbolToAtomicNumber(material)
        density = xraylib.ElementDensity(Z)
        delta = 1.0-xraylib.Refractive_Index_Re(material, energy, density)
        
        # Add the delta value to the delta list
        delta_list.append(delta)
    
    return np.array(delta_list)
    
def materials_to_linear_attenuation(material_list, energy):
    """
    Description:
        Convert a list of material names to a list of linear attenuation 
        coefficients at a given energy.
    
    Parameters:
        material_list (list): A list of material names.
        energy (float): The energy in keV.
    
    Returns:
        list: A list of linear attenuation coefficient values (in m^-1) for the given materials at the given energy.
    """

    # The list to store linear attenuation coefficient values
    mu_list = []

    # Iterate through each material in the input list
    for material in material_list:
        # Compute the delta value for the current material at the given energy
        Z = xraylib.SymbolToAtomicNumber(material)
        density = xraylib.ElementDensity(Z)
        # Compute the mass attenuation coefficient in cm^2/g
        #mass_attenuation = xraylib.CS_Photo(Z, energy)
        mass_attenuation = xraylib.CS_Total(Z, energy)
        # Convert mass attenuation coefficient to linear attenuation coefficient in m^-1
        mu = mass_attenuation * density * 100.0       
        # Add the linear attenuation coefficient value to the list
        mu_list.append(mu)
    
    return np.asarray(mu_list)

def absorptionaperture(x, n1mud, sigma, n1mur):
    '''
    Description:
        Calculates absorption aperture (See XS for more)
    
    Parameters:
        x       :
        n1mud   :
        sigma   :
        n1mur   :
        
    Returns:
        absorption aperture?
        
    '''
    
    numerator = np.exp(-(x**2/(2*sigma**2))) * np.exp(-n1mur*(x**2) - n1mud)
    denominator = np.exp(-n1mud)
    return numerator / denominator - 0.5


def find_levels(array, levels, direction='forward'):
    """
    Find the first indices at which the array crosses specified levels and the corresponding crossed values.

    Parameters:
        array (numpy.ndarray): An array of numbers.
        levels (float or numpy.ndarray): A number or an array of levels to find crossings.
        direction (str, optional): The searching direction. Defaults to 'forward'.
                                   Can be either 'forward' or 'backward'.

    Returns:
        tuple: A tuple containing two arrays:
            - An array of first indices at which the array crosses the specified levels.
            - An array of first crossed values at the corresponding indices.
    """

    # Convert a single level to a numpy array
    if isinstance(levels, (int, float)):
        levels = np.array([levels])

    indices = []
    values = []

    # Compute the max and min of the array ignoring NaNs
    max_val = np.nanmax(array)
    min_val = np.nanmin(array)

    for level in levels:
        # If level is out of bounds
        if level > max_val or level < min_val:
            indices.append(-1)
            values.append(np.nan)
            continue

        crossings = []

        if direction == 'forward':
            for i in range(1, len(array)):
                if np.isnan(array[i - 1]) or np.isnan(array[i]):
                    continue
                if (array[i - 1] < level <= array[i]) or (array[i - 1] > level >= array[i]):
                    crossings.append(i - 1)
                    break
        elif direction == 'backward':
            for i in range(len(array) - 2, -1, -1):
                if np.isnan(array[i + 1]) or np.isnan(array[i]):
                    continue
                if (array[i + 1] < level <= array[i]) or (array[i + 1] > level >= array[i]):
                    crossings.append(i)
                    break
        else:
            raise ValueError("Invalid direction. It should be either 'forward' or 'backward'.")

        if len(crossings) > 0:
            idx = crossings[0]
            indices.append(idx)
            values.append(array[idx])
        else:
            # In case no crossing is found within the range
            indices.append(-1)
            values.append(np.nan)

    return np.array(indices), np.array(values)

def calc_tf1_data(num_configs, radii, materials, energy_keV, wl, numlens, 
                  lens_loc, beam, bl, crl, slit_H, slit_V, thickerr, 
                  flag_HE = False, verbose = False):

    '''
    Description:
        Begains base calculation for first tranfocator lookup table.  Returns
        quantities needed by three system types (single transfocator, double
        transfocator, and transfocator + KB mirror
    Parameters:
        num_configs : integer
            number of configs for system (i.e. entries in lookup table)
        radii       : list of float
            radius of each lens in transfocator (in m)
        materials   : list of strings
            lens material for each stack
        energy_keV  : float
            beam energy in keV
        wl          : float
            wavelength of beam photon (in m)
        numlens     : list of int
            number of lenses in each stack
        lens_loc    : list of float
            position of each stack with respect to transfocator center
        beam        : dictionary
            beam properties
        bl          : dictionary
            beamline properties (element locations)
        crl         : dictionary
            Transfocator properties
        slit_H      : float
            horizontal slit size
        slit_V      : float
            vertical slit size
        thickerr:   : list of float
            thickness error for each stack
        flag_HE     : boolean 
            flag on whether to use thickness error in focal size calc (True) or 
            not (False)
        verbose     : boolean               
            flag for printing to iocConsole
            
    Returns:
        Dictionary with the following keys
            q1_list                   : numpy array of float
                image position for each configuration relative to element center
            dq1_list                  : numpy array of float
                image position for each configuration measured from source
            aperL1H_list              : numpy array of float
                absorption H aperture of transfocator 1 for all configurations
            aperL1V_list              : numpy array of float
                absorption V aperture of transfocator 1 for all configurations
            FWHM1H_list               : numpy array of float
                H focal size at the focus of transfocator 1
            FWHM1V_list               : numpy array of float
                V focal size at the focus of transfocator 1
            L1_invF_list_sorted       : numpy array of float
                equivalent 1/f for each configuration (sorted by increasing value)
            L1_invF_list_sort_indices : numpy array of float
                sorted indices for L1_invF_list_sorted
                
    '''


    if verbose:
        print(30*'*')
        print(f'Energy: {energy_keV} keV')
        print(f'Hor slit size: {slit_H} m')
        print(f'Ver slit size: {slit_V} m')
        print(30*'*')

          
    lookupTable = np.empty(num_configs)
            
    sigmaH = beam['sigmaH']
    sigmaV = beam['sigmaV']
    sigmaHp = beam['sigmaHp']
    sigmaVp = beam['sigmaVp']

    d_StoL1 = bl['d_StoL1']
    d_Stof = bl['d_Stof']

    d_min = crl['d_min']

   
    L1_D        = np.asarray([lookup_diameter(rad) for rad in radii])    # CRL1 diameters for each stack
    L1_delta    = materials_to_deltas(materials, energy_keV)             # delta values for CRL1 stacks
    L1_mu       = materials_to_linear_attenuation(materials, energy_keV) # mu values for CRL1 stacks
    L1_Feq      = radii/(2*np.asarray(numlens)*L1_delta) + lens_loc                       # CRL1 equivalent f in m for each stack
    L1_index_n  = 2**L1_Feq.size                                        # Total number of combinations for CRL1
    L1_invF_list= np.zeros(L1_index_n)                                  

    # List of equivalent 1/f in m^-1 for CRL1
    L1_invF_list = np.asarray([sum(index_to_binary_list(i, L1_Feq.size)/L1_Feq) for i in range(L1_index_n)]) 
    # Sort the L1_invF list inverse of focal length
    L1_invF_list_sort_indices = np.argsort(L1_invF_list)
    L1_invF_list_sorted       = L1_invF_list[L1_invF_list_sort_indices]
    # image position of CRL1 for all configurations (sorted by inverse focal length)
    q1_list  = 1/(L1_invF_list_sorted - 1/d_StoL1)      
    dq1_list = q1_list - (d_Stof - d_StoL1)
        
    # Start generating focal size list as a function of CRL1 setting
    sigma1H         = (sigmaH**2 + (sigmaHp*d_StoL1)**2)**0.5   # sigma beam size before CRL1
    sigma1V         = (sigmaV**2 + (sigmaVp*d_StoL1)**2)**0.5   # sigma beam size before CRL1
    L1_n1mud_list   = np.zeros(L1_index_n)                      # List of n1*mu*d_min for all possible CRL1 configurations
    L1_n1muR_list   = np.zeros(L1_index_n)                      # List of n1*mu/R for all possible CRL1 configurations
    aperL1H_list    = np.zeros(L1_index_n)                      # absorption H aperture of CRL1 for all configurations
    aperL1V_list    = np.zeros(L1_index_n)                      # absorption V aperture of CRL1 for all configurations
    diameter1_list  = np.zeros(L1_index_n)                      # CRL1 diameter for all possible configurations
    FWHM1H_list     = np.zeros(L1_index_n)                      # H focal size at the focus of CRL1
    FWHM1V_list     = np.zeros(L1_index_n)                      # V focal size at the focus of CRL1
    Strehl_list     = np.zeros(L1_index_n)                      # Strehl ratio based on lens thickness error

    for i in range(L1_index_n):
        # absorption aperture is a function of CRL absorption/physical aperture, incident beam size, and physical slits
        L1_n1mud_list[i] = np.sum(index_to_binary_list(L1_invF_list_sort_indices[i], L1_Feq.size)*np.array(L1_mu*numlens*d_min))
#        L1_n1mud_list[i] = np.sum(index_to_binary_list(i, L1_Feq.size)*np.array(L1_mu*numlens*d_min))
        L1_n1muR_list[i] = np.sum(index_to_binary_list(L1_invF_list_sort_indices[i], L1_Feq.size)*np.array(L1_mu*numlens/radii))
#        L1_n1muR_list[i] = np.sum(index_to_binary_list(i, L1_Feq.size)*np.array(L1_mu*numlens/radii))
        solution = root_scalar(absorptionaperture, args=(L1_n1mud_list[i], sigma1H, L1_n1muR_list[i]), bracket=[0.0, 2*sigma1H], method='bisect')
        aperL1H_list[i] = solution.root*2.0
        solution = root_scalar(absorptionaperture, args=(L1_n1mud_list[i], sigma1V, L1_n1muR_list[i]), bracket=[0.0, 2*sigma1V], method='bisect')
        aperL1V_list[i] = solution.root*2.0
        mask = (np.array(index_to_binary_list(L1_invF_list_sort_indices[i], L1_Feq.size)) == 1)
#        mask = (np.array(index_to_binary_list(i, L1_Feq.size)) == 1)
        if np.all(mask == False):
            diameter1_list[i] = np.inf
        else:
            diameter1_list[i] = np.min(L1_D[mask])
        aperL1H_list[i] = min(aperL1H_list[i], diameter1_list[i], slit_H)
        aperL1V_list[i] = min(aperL1V_list[i], diameter1_list[i], slit_V)
        phase_error_tmp = np.linalg.norm(index_to_binary_list(L1_invF_list_sort_indices[i], L1_Feq.size)*np.array(thickerr*L1_delta)*2*np.pi/wl)
#        phase_error_tmp = np.linalg.norm(index_to_binary_list(i, L1_Feq.size)*np.array(thickerr*L1_delta)*2*np.pi/wl)
        Strehl_list[i] = np.exp(-phase_error_tmp**2)

        
    # FWHMbeam size at CRL1 focus
    FWHM1H_list  = ((0.88*wl*q1_list/aperL1H_list)**2 + (2.355*sigmaH*q1_list/d_StoL1)**2)**0.5
    FWHM1V_list  = ((0.88*wl*q1_list/aperL1V_list)**2 + (2.355*sigmaV*q1_list/d_StoL1)**2)**0.5
    if flag_HE:
        FWHM1H_list *= (Strehl_list)**(-0.5)
        FWHM1V_list *= (Strehl_list)**(-0.5)    

    return {'q1_list': q1_list, 'dq1_list': dq1_list, 'aperL1H_list': aperL1H_list, 
            'aperL1V_list': aperL1V_list, 'FWHM1H_list': FWHM1H_list, 'FWHM1V_list': FWHM1V_list,
            'L1_invF_list_sort_indices': L1_invF_list_sort_indices, 
            'L1_invF_list_sorted': L1_invF_list_sorted, 'L1_index_n': L1_index_n}


def calc_1x_lu_table(num_configs, radii, materials, energy_keV, wl, numlens, 
                     lens_loc, beam, bl, crl, slit_H, slit_V, thickerr, 
                     flag_HE = False, verbose = False):
    '''
    Description:
        Lookup table calculation for single CRL system
        
    Parameters:
        num_configs     : number of CRL1 configurations
        radii           : List of lens radii in CRL1
        materials       : List of lens materials in CRL1
        energy_keV      : incident beam energy in keV
        wl              : beam wavelength in nm
        numlens         : number of lenses in CRL 1 (list)
        lens_loc        : lens locations wrt to CRL 1 center 
        beam            : beam properties dictionary
        bl              : beamline properties dictionary
        crl             : CRL properties dictionary
        slit_H          : float
            horizontal slit size
        slit_V          : float
            vertical slit size
        thickerr        : thickness error
        flag_HE         : Flag to include thickness error in calculation
        verbose         : Flag to print messages to IOC console
    Returns:
        FWHM_atsample           : focal size in meters
        invF_list_sort_indices  : elements are n-bit config for CRL1, sorted by increasing equivalent 1/f
        invF_list_sorted        : List of equivalent 1/f in m^-1 for CRL1, sorted by increasing value
    '''    
    data_dict = calc_tf1_data(num_configs, radii, materials, energy_keV, wl, numlens, 
                     lens_loc, beam, bl, crl['1'], slit_H, slit_V, thickerr, 
                     flag_HE = flag_HE, verbose = verbose)

    q1_list = data_dict['q1_list']
    dq1_list = data_dict['dq1_list']
    aperL1H_list = data_dict['aperL1H_list']
    aperL1V_list = data_dict['aperL1V_list']
    FWHM1H_list = data_dict['FWHM1H_list']
    FWHM1V_list = data_dict['FWHM1V_list']
    L1_invF_list_sort_indices = data_dict['L1_invF_list_sort_indices']
    L1_invF_list_sorted = data_dict['L1_invF_list_sorted']

    FWHM1H_atsample_list = (FWHM1H_list**2 + (aperL1H_list*dq1_list/q1_list)**2)**0.5
    FWHM1V_atsample_list = (FWHM1V_list**2 + (aperL1V_list*dq1_list/q1_list)**2)**0.5
    FWHM_atsample_list   = (FWHM1H_atsample_list*FWHM1V_atsample_list)**0.5

    invF_list_sort_indices = {'1': L1_invF_list_sort_indices, '2': None}
    invF_list_sorted = {'1': L1_invF_list_sorted, '2': None}

    return FWHM_atsample_list, invF_list_sort_indices, invF_list_sorted


def calc_2x_lu_table(num_configs, L1_radii, L1_materials, L2_radii, L2_materials, 
                     energy_keV, wl, numlens, L1_lens_loc, L2_lens_loc, beam, bl, 
                     crl, slits, L1_thickerr, L2_thickerr, 
                     flag_HE = False, verbose = False):

    '''
    Description:
        Lookup table calculation for Double CRL system
        
    Parameters:
        num_configs     : number of CRL1 configurations
        L1_radii        : List of lens radii in CRL1
        L1_materials    : List of lens materials in CRL1
        L2_radii        : List of lens radii in CRL2
        L2_materials    : List of lens materials in CRL2
        energy_keV      : incident beam energy in keV
        wl              : beam wavelength in nm
        numlens         : number of lenses in each crl dictionary of lists
        L1_lens_loc     : lens locations wrt to CRL 1 center 
        L2_lens_loc     : lens locations wrt to CRL 1 center
        beam            : beam properties dictionary
        bl              : beamline properties dictionary
        crl             : CRL properties dictionary of dictionaries
        slits           : slits sizes dictionary
        L1_thickerr     : CRL1 thickness errors
        L2_thickerr     : CRL2 thickness errors
        flag_HE         : Flag to include thickness error in calculation
        verbose         : Flag to print messages to IOC console
    Returns:
        FWHM_atsample           : focal size in meters
        invF_list_sort_indices  : dictionary (L1, L2) for n-bit configs for each 
                                  CRL, sorted by increasing equivalent 1/f
        invF_list_sorted        : dictionary (L1, L2) for equivalent 1/f in m^-1 
                                  for each CRL, sorted by increasing value
        invf2_indices           : each element is the CRL2 n-bit configuration 
                                  and each index is the sorted index for CRL1, 
                                  e.g. 
    '''


    d_StoL2 = bl['d_StoL2']
    d_StoL1 = bl['d_StoL1']
    d_Stof = bl['d_Stof']

    data_dict = calc_tf1_data(num_configs, L1_radii, L1_materials, energy_keV, wl, numlens['1'], 
                     L1_lens_loc, beam, bl, crl['1'], slits['1']['hor'], slits['1']['vert'], 
                     L1_thickerr, 
                     flag_HE = flag_HE, verbose = verbose)

    q1_list = data_dict['q1_list']
    dq1_list = data_dict['dq1_list']
    aperL1H_list = data_dict['aperL1H_list']
    aperL1V_list = data_dict['aperL1V_list']
    FWHM1H_list = data_dict['FWHM1H_list']
    FWHM1V_list = data_dict['FWHM1V_list']
    L1_invF_list_sort_indices = data_dict['L1_invF_list_sort_indices']
    L1_invF_list_sorted = data_dict['L1_invF_list_sorted']
    L1_index_n = data_dict['L1_index_n']

    L2_D        = np.asarray([lookup_diameter(rad) for rad in L2_radii])# CRL2 diameters for each stack
    L2_delta    = materials_to_deltas(L2_materials, energy_keV)             # Delta values for CRL2 stacks
    L2_mu       = materials_to_linear_attenuation(L2_materials, energy_keV) # mu values for CRL2 stacks
    L2_Feq      = L2_radii/(2*np.asarray(numlens['2'])*L2_delta)+L2_lens_loc                         # CRL2 equivalent f in m for each stack
    L2_index_n   = 2**L2_Feq.size                               # Total number of combinations for CRL2


    # List of equivalent 1/f in m^-1 for CRL2
    L2_invF_list = np.asarray([sum(index_to_binary_list(i, L2_Feq.size)/L2_Feq) for i in range(L2_index_n)])

    # Sort the L2_invF list (to avoid zigzagging)
    L2_invF_list_sort_indices = np.argsort(L2_invF_list)
    L2_invF_list_sorted       = L2_invF_list[L2_invF_list_sort_indices]

#    sigma2H_list    = np.zeros(L1_index_n)                      # sigma beam size before CRL2
#    sigma2V_list    = np.zeros(L1_index_n)                      # sigma beam size before CRL2

    # Sigma beam size before CRL2
    sigma2H_list = (((0.88*wl*(d_StoL2-d_StoL1))/aperL1H_list)**2 + (aperL1H_list*(1-(d_StoL2-d_StoL1)/q1_list))**2)**0.5/2.355
    sigma2V_list = (((0.88*wl*(d_StoL2-d_StoL1))/aperL1V_list)**2 + (aperL1V_list*(1-(d_StoL2-d_StoL1)/q1_list))**2)**0.5/2.355

    p2_list      = d_StoL2 - d_StoL1 - q1_list           # p2 for CRL2 for all possible CRL1 configurations
    invf2_list   = 1.0/p2_list + 1/(d_Stof - d_StoL2)    # f2^-1 for CRL2 to match CRL1 for all possible CRL1 configurations
    #L2_config_index   = np.zeros(L1_index_n)            # CRL2 configueration index to match CRL1

    #invf2_indices, invf2_values = find_closest_values_auto(L2_invF_list_sorted, invf2_list)
    #invf2_indices, invf2_values = find_levels_left(L2_invF_list_sorted, invf2_list)
    invf2_indices, invf2_values = find_levels(L2_invF_list_sorted, invf2_list, direction = 'forward')

    nan_positions = np.where(invf2_indices == -1)
    invf2_values[nan_positions] = np.nan                # only f2^-1 values that can be matched with CRL1
    q2_list  = 1/(invf2_values - 1/p2_list)
    dq2_list = q2_list - (d_Stof - d_StoL2)

    L2_n2mud_list   = np.zeros(L1_index_n)              # List of n2*mu*d_min for all possible CRL1 configurations
    L2_n2muR_list   = np.zeros(L1_index_n)              # List of n2*mu/R for all possible CRL1 configurations
    aperL2H_list    = np.zeros(L1_index_n)              # absorption H aperture of CRL2 for all CRL1 configurations
    aperL2V_list    = np.zeros(L1_index_n)              # absorption V aperture of CRL2 for all CRL1 configurations
    diameter2_list  = np.zeros(L1_index_n)              # CRL2 diameter for all possible CRL1 configurations
    FWHM2H_list     = np.zeros(L1_index_n)              # H focal size at the focus of CRL2 matching all possible CRL1 configurations
    FWHM2V_list     = np.zeros(L1_index_n)              # V focal size at the focus of CRL2 matching all possible CRL1 configurations
    FWHM_list       = np.zeros(L1_index_n)              # Focal size sqrt(H*V) at the focus of CRL2 matching all possible CRL1 configurations
    Strehl2_list     = np.zeros(L1_index_n)                      # Strehl ratio based on lens thickness error

    for i in range(L1_index_n):
        if invf2_indices[i] != -1:
            # absorption aperture is a function of CRL absorption/physical aperture, incident beam size, and physical slits
            L2_n2mud_list[i] = np.sum(index_to_binary_list(L2_invF_list_sort_indices[invf2_indices[i]], L2_Feq.size)*np.array(L2_mu*numlens['2']*crl['2']['d_min']))
            L2_n2muR_list[i] = np.sum(index_to_binary_list(L2_invF_list_sort_indices[invf2_indices[i]], L2_Feq.size)*np.array(L2_mu*numlens['2']/L2_radii))
            solution = root_scalar(absorptionaperture, args=(L2_n2mud_list[i], sigma2H_list[i], L2_n2muR_list[i]), bracket=[0.0, 2*sigma2H_list[i]], method='bisect')
            aperL2H_list[i] = solution.root*2.0
            solution = root_scalar(absorptionaperture, args=(L2_n2mud_list[i], sigma2V_list[i], L2_n2muR_list[i]), bracket=[0.0, 2*sigma2V_list[i]], method='bisect')
            aperL2V_list[i] = solution.root*2.0
            mask = (np.array(index_to_binary_list(L2_invF_list_sort_indices[invf2_indices[i]], L2_Feq.size)) == 1)
            if np.all(mask == False):
                diameter2_list[i] = np.inf
            else:
                diameter2_list[i] = np.min(L2_D[mask])
            aperL2H_list[i] = min(aperL2H_list[i], diameter2_list[i], slits['2']['hor'])
            aperL2V_list[i] = min(aperL2V_list[i], diameter2_list[i], slits['2']['vert'])
            phase_error_tmp = np.linalg.norm(index_to_binary_list(L2_invF_list_sort_indices[invf2_indices[i]], L2_Feq.size)*np.array(L2_thickerr*L2_delta)*2*np.pi/wl)
            Strehl2_list[i] = np.exp(-phase_error_tmp**2)
    aperL2H_list[nan_positions] = np.nan
    aperL2V_list[nan_positions] = np.nan
    Strehl2_list[nan_positions] = np.nan

    # FWHMbeam size at focus
    FWHM2H_list = ((0.88*wl*q2_list/aperL2H_list)**2 + (FWHM1H_list*q2_list/p2_list)**2)**0.5
    FWHM2V_list = ((0.88*wl*q2_list/aperL2V_list)**2 + (FWHM1V_list*q2_list/p2_list)**2)**0.5
    if flag_HE:
        FWHM2H_list *= (Strehl2_list)**(-0.5)
        FWHM2V_list *= (Strehl2_list)**(-0.5)
    FWHM_list   = (FWHM2H_list*FWHM2V_list)**0.5

    # FWHMbeam size at sample
    FWHM2H_atsample_list = (FWHM2H_list**2 + (aperL2H_list*dq2_list/q2_list)**2)**0.5
    FWHM2V_atsample_list = (FWHM2V_list**2 + (aperL2V_list*dq2_list/q2_list)**2)**0.5
    FWHM_atsample_list   = (FWHM2H_atsample_list*FWHM2V_atsample_list)**0.5
    
    invF_list_sort_indices = {'1': L1_invF_list_sort_indices, '2': L2_invF_list_sort_indices}
    invF_list_sorted = {'1': L1_invF_list_sorted, '2': L2_invF_list_sorted}

    return FWHM_atsample_list, invF_list_sort_indices, invF_list_sorted, invf2_indices



def calc_kb_lu_table(num_configs, radii, materials, energy_keV, wl, numlens, 
                      lens_loc, beam, bl, crl, kb, slits, thickerr, 
                      flag_HE = False, verbose = False):
    '''
    Description:
        Lookup table calculation for CRL + KB system
        
    Parameters:
        num_configs     : number of CRL1 configurations
        radii           : List of lens radii in CRL1
        materials       : List of lens materials in CRL1
        energy_keV      : incident beam energy in keV
        wl              : beam wavelength in nm
        numlens         : number of lenses in CRL 1 (list)
        lens_loc        : lens locations wrt to CRL 1 center 
        beam            : beam properties dictionary
        bl              : beamline properties dictionary
        crl             : CRL properties dictionary
        kb              : KB properties dictionary
        slits           : slits sizes dictionary
        thickerr        : thickness error
        flag_HE         : Flag to include thickness error in calculation
        verbose         : Flag to print messages to IOC console
    Returns:
        FWHM_atsample           : focal size in meters
        invF_list_sort_indices  :
        invF_list_sorted        :
    '''

    d_StoL1 = bl['d_StoL1']
    d_Stof = bl['d_Stof']

    KBH_q = kb['KBH_q']
    KBH_L = kb['KBH_L']
    KBH_slit = slits['kb']['hor']
    KBH_p_limit = kb['KBH_p_limit']

    KBV_q = kb['KBV_q']
    KBV_L = kb['KBV_L']
    KBV_slit = slits['kb']['vert']
    KBV_p_limit = kb['KBV_p_limit']

    KB_theta = kb['KB_theta']

    data_dict = calc_tf1_data(num_configs, radii, materials, energy_keV, wl, numlens['1'], 
                     lens_loc, beam, bl, crl, slits['1']['hor'], slits['1']['vert'], thickerr, 
                     flag_HE = flag_HE, verbose = verbose)

    q1_list = data_dict['q1_list']
    dq1_list = data_dict['dq1_list']
    aperL1H_list = data_dict['aperL1H_list']
    aperL1V_list = data_dict['aperL1V_list']
    FWHM1H_list = data_dict['FWHM1H_list']
    FWHM1V_list = data_dict['FWHM1V_list']
    L1_invF_list_sort_indices = data_dict['L1_invF_list_sort_indices']
    L1_invF_list_sorted = data_dict['L1_invF_list_sorted']

    aperL2H_list    = np.zeros(L1_index_n)              # absorption H aperture of KBH for all CRL1 configurations
    aperL2V_list    = np.zeros(L1_index_n)              # absorption V aperture of KBV for all CRL1 configurations
    FWHM2H_list     = np.zeros(L1_index_n)              # H focal size at the focus of KBH matching all possible CRL1 configurations
    FWHM2V_list     = np.zeros(L1_index_n)              # V focal size at the focus of KBV matching all possible CRL1 configurations


    # Sigma beam size before KB
    sigma2H_list  = (((0.88*wl*(d_Stof - KBH_q - d_StoL1))/aperL1H_list)**2 + (aperL1H_list*(1-(d_Stof-KBH_q-d_StoL1)/q1_list))**2)**0.5/2.355
    sigma2V_list  = (((0.88*wl*(d_Stof - KBV_q - d_StoL1))/aperL1V_list)**2 + (aperL1V_list*(1-(d_Stof-KBV_q-d_StoL1)/q1_list))**2)**0.5/2.355

    KBH_p_list     = d_Stof - KBH_q - d_StoL1 - q1_list           # p for KBH for all possible CRL1 configurations
    KBV_p_list     = d_Stof - KBV_q - d_StoL1 - q1_list           # p for KBH for all possible CRL1 configurations
    nan_positionsH = np.where((KBH_p_list > -KBH_p_limit) & (KBH_p_list < KBH_p_limit))   # p too close to mirror
    KBH_p_list[nan_positionsH] = np.nan
    nan_positionsV = np.where((KBV_p_list > -KBV_p_limit) & (KBV_p_list < KBV_p_limit))
    KBV_p_list[nan_positionsV] = np.nan

    KBH_R_list = 2/np.sin(KB_theta)/(1/KBH_p_list+1/KBH_q)
    KBV_R_list = 2/np.sin(KB_theta)/(1/KBV_p_list+1/KBV_q)
    aperL2H_list = np.minimum(sigma2H_list*2.355, np.minimum(KBH_L*KB_theta, KBH_slit))
    aperL2V_list = np.minimum(sigma2V_list*2.355, np.minimum(KBV_L*KB_theta, KBV_slit))

    # FWHMbeam size at focus (coincident with sample for KB case)
    FWHM2H_list = ((0.88*wl*KBH_q/aperL2H_list)**2 + (FWHM1H_list*KBH_q/KBH_p_list)**2)**0.5
    FWHM2V_list = ((0.88*wl*KBV_q/aperL2V_list)**2 + (FWHM1V_list*KBV_q/KBV_p_list)**2)**0.5
    FWHM_list   = (FWHM2H_list*FWHM2V_list)**0.5

    invF_list_sort_indices = {'1': L1_invF_list_sort_indices, '2': None}
    invF_list_sorted = {'1': L1_invF_list_sorted, '2': None}

    return FWHM_atsample_list, invF_list_sort_indices, invF_list_sorted


def calc_2xCRL_focus(index1, index2, L1_radii, L1_materials, L2_radii, L2_materials, energy_keV, wl, numlens, 
                      L1_lens_loc, L2_lens_loc, beam, bl, crl, slits, L1_thickerr, L2_thickerr,
                      flag_HE = False, verbose = False):
    '''
    Description:
        When 2nd CRL is "tweaked", system configuration is no longer on lookup 
        table and need to calculate focus for new configuration (index1, index2)
        
    Parameters:
        index1          : TF1 n-bit configuration
        index2          : TF2 n-bit configuration
        L1_radii        : List of lens radii in CRL1
        L1_materials    : List of lens materials in CRL1
        L2_radii        : List of lens radii in CRL2
        L2_materials    : List of lens materials in CRL2
        energy_keV      : incident beam energy in keV
        wl              : beam wavelength in nm
        numlens         : number of lenses in each crl (list)
        L1_lens_loc     : lens locations wrt to CRL 1 center 
        L2_lens_loc     : lens locations wrt to CRL 1 center
        beam            : beam properties dictionary
        bl              : beamline properties dictionary
        crl             : CRL properties list of dictionaries
        slits           : slits sizes dictionary
        L1_thickerr     : CRL1 thickness errors
        L2_thickerr     : CRL2 thickness errors
        flag_HE         : Flag to include thickness error in calculation
        verbose         : Flag to print messages to IOC console
    Returns:
        FWHM_atsample   : focal size in meters
    '''

    sigmaH = beam['sigmaH']
    sigmaV = beam['sigmaV']
    sigmaHp = beam['sigmaHp']
    sigmaVp = beam['sigmaVp']

    d_StoL2 = bl['d_StoL2']
    d_StoL1 = bl['d_StoL1']
    d_Stof = bl['d_Stof']

        
    L1_D        = np.asarray([lookup_diameter(rad) for rad in L1_radii])    # CRL1 diameters for each stack
    L1_delta    = materials_to_deltas(L1_materials, energy_keV)             # delta values for CRL1 stacks
    L1_mu       = materials_to_linear_attenuation(L1_materials, energy_keV) # mu values for CRL1 stacks
    L1_Feq      = L1_radii/(2*np.asarray(numlens['1'])*L1_delta) + L1_lens_loc                       # CRL1 equivalent f in m for each stack

    L2_D        = np.asarray([lookup_diameter(rad) for rad in L2_radii])# CRL2 diameters for each stack
    L2_delta    = materials_to_deltas(L2_materials, energy_keV)             # Delta values for CRL2 stacks
    L2_mu       = materials_to_linear_attenuation(L2_materials, energy_keV) # mu values for CRL2 stacks
    L2_Feq      = L2_radii/(2*np.asarray(numlens['2'])*L2_delta) + L2_lens_loc                       # CRL2 equivalent f in m for each stack

    # Calculation block
    L1_invF = np.sum(index_to_binary_list(index1, L1_Feq.size)/L1_Feq)  # f^-1 for CRL1
    L2_invF = np.sum(index_to_binary_list(index2, L2_Feq.size)/L2_Feq)  # f^-1 for CRL2
    q1      = 1/(L1_invF - 1/d_StoL1)                                   # focal position of CRL1
    sigma1H = (sigmaH**2 + (sigmaHp*d_StoL1)**2)**0.5                   # sigma beam size before CRL1
    sigma1V = (sigmaV**2 + (sigmaVp*d_StoL1)**2)**0.5                   # sigma beam size before CRL1

    # absorption aperture is a function of CRL absorption/physical aperture, incident beam size, and physical slits
    L1_n1mud = np.sum(index_to_binary_list(index1, L1_Feq.size)*np.array(L1_mu*numlens['1']*crl['1']['d_min']))
    L1_n1muR = np.sum(index_to_binary_list(index1, L1_Feq.size)*np.array(L1_mu*numlens['1']/L1_radii))
    solution = root_scalar(absorptionaperture, args=(L1_n1mud, sigma1H, L1_n1muR), bracket=[0.0, 2*sigma1H], method='bisect')
    aperL1H  = solution.root*2.0
    solution = root_scalar(absorptionaperture, args=(L1_n1mud, sigma1V, L1_n1muR), bracket=[0.0, 2*sigma1V], method='bisect')
    aperL1V  = solution.root*2.0
    mask     = (np.array(index_to_binary_list(index1, L1_Feq.size)) == 1)
    if np.all(mask == False):
        diameter1 = np.inf
    else:
        diameter1 = np.min(L1_D[mask])
    aperL1H = min(aperL1H, diameter1, slits['1']['hor'])
    aperL1V = min(aperL1V, diameter1, slits['1']['vert'])
    phase_error_tmp1 = np.linalg.norm(index_to_binary_list(index1, L1_Feq.size)*np.array(L1_thickerr*L1_delta)*2*np.pi/wl)
    Strehl1 = np.exp(-phase_error_tmp1**2)

    # FWHMbeam size at CRL1 focus
    FWHM1H  = ((0.88*wl*q1/aperL1H)**2 + (2.355*sigmaH*q1/d_StoL1)**2)**0.5
    FWHM1V  = ((0.88*wl*q1/aperL1V)**2 + (2.355*sigmaV*q1/d_StoL1)**2)**0.5
    if flag_HE:
        FWHM1H *= (Strehl1)**(-0.5)
        FWHM1V *= (Strehl1)**(-0.5)

    # Sigma beam size before CRL2
    sigma2H = (((0.88*wl*(d_StoL2-d_StoL1))/aperL1H)**2 + (aperL1H*(1-(d_StoL2-d_StoL1)/q1))**2)**0.5/2.355
    sigma2V = (((0.88*wl*(d_StoL2-d_StoL1))/aperL1V)**2 + (aperL1V*(1-(d_StoL2-d_StoL1)/q1))**2)**0.5/2.355

    p2      = d_StoL2 - d_StoL1 - q1    # p2 for CRL2
    q2      = 1/(L2_invF - 1/p2)        # q2 for CRL2 calculated from CRL2 index and p2
    dq2     = q2 - (d_Stof - d_StoL2)   # off focus distance

    L2_n2mud = np.sum(index_to_binary_list(index2, L2_Feq.size)*np.array(L2_mu*numlens['2']*crl['2']['d_min']))
    L2_n2muR = np.sum(index_to_binary_list(index2, L2_Feq.size)*np.array(L2_mu*numlens['2']/L2_radii))
    solution = root_scalar(absorptionaperture, args=(L2_n2mud, sigma2H, L2_n2muR), bracket=[0.0, 2*sigma2H], method='bisect')

    aperL2H = solution.root*2.0
    solution = root_scalar(absorptionaperture, args=(L2_n2mud, sigma2V, L2_n2muR), bracket=[0.0, 2*sigma2V], method='bisect')
    aperL2V = solution.root*2.0
    mask = (np.array(index_to_binary_list(index2, L2_Feq.size)) == 1)
    if np.all(mask == False):
        diameter2 = np.inf
    else:
        diameter2 = np.min(L2_D[mask])
    aperL2H = min(aperL2H, diameter2, slits['2']['hor'])
    aperL2V = min(aperL2V, diameter2, slits['2']['vert'])
    phase_error_tmp2 = np.linalg.norm(index_to_binary_list(index2, L2_Feq.size)*np.array(L2_thickerr*L2_delta)*2*np.pi/wl)
    Strehl2 = np.exp(-phase_error_tmp2**2)

    FWHM2H = ((0.88*wl*q2/aperL2H)**2 + (FWHM1H*q2/p2)**2)**0.5
    FWHM2V = ((0.88*wl*q2/aperL2V)**2 + (FWHM1V*q2/p2)**2)**0.5
    if flag_HE:
        FWHM2H *= (Strehl2)**(-0.5)
        FWHM2V *= (Strehl2)**(-0.5)

    FWHM2H_atsample = (FWHM2H**2 + (aperL2H*dq2/q2)**2)**0.5
    FWHM2V_atsample = (FWHM2V**2 + (aperL2V*dq2/q2)**2)**0.5
    FWHM_atsample   = (FWHM2H_atsample*FWHM2V_atsample)**0.5

    return FWHM_atsample
    
