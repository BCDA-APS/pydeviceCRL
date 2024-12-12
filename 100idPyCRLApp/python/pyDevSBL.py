import numpy as np
import os

class surrogateBL():
    
    def __init__(self, beamline):
        '''
        beamline : Shadow-based beamline
        '''
           
        #create necessary inputs/output; initialize to zero/empty
        print(80*'#')
        print('Initializing pyDevice surrogate beamline object')
        self.beamline = beamline
        print(80*'#')
        self.motorPositions = []
        self.verbose = 1   

    def setupMotorLimits(self, subsFile):
        '''
        Go through shadow beamline and get DOF limits and set them as the motor
        limits
        
        Don't have a good way to do this; will manually update simulated motor
        substitutions file with limits
        '''
        pass
        
       
    def elementUpdate(self, element, dof, position):
        '''
        Update element position array but do not trigger detector
        '''
        if self.verbose: print(f'Setting element {element}, DOF {dof} to {position}')
        self.beamline.pos[element, dof]=position
        if self.verbose: print(f'Setting complete')
     
        
    def detectorUpdate(self, zeroD = True, nbins = 64):
        '''
        Update detector(s)
        '''
        
        # Update optical element positions
        if self.verbose: print(f'Adjusting OE settings ')
        self.beamline.adjust() 
        if self.verbose: print(f'Adjusting OE settings complete')
        
        # Trigger ray tracing
        if self.verbose: print(f'Running Shadow simulation')        
        if zeroD:
            results = self.beamline.run()
        else:
            results = self.beamline.run(nbins = nbins, nolost=1) 
        if self.verbose: 
            print(f'Run complete')
            print(f'Results: {results}')

        if zeroD:
            self.intensity = int(results)       
            pydev.iointr('new_intensity', self.intensity)
        else:     
            self.image_data = results['histogram']
            self.flattened_data = self.image_data.flatten() # need to flatten for waveform
            self.h_edges = results['bin_h_edges']
            self.v_edges = results['bin_v_edges']
            pydev.iointr('raw_image', self.flattened_data.tolist())
            pydev.iointr('h_coord', self.h_edges.tolist())
            pydev.iointr('v_coord', self.v_edges.tolist())
            self.intensity = int(results['good_rays'])
            pydev.iointr('new_intensity', self.intensity)
            
    def updateVerbosity(self, verbosity):
        '''
        Turn on minor printing
        '''
        print(f'Verbosity set to {verbosity}')
        self.verbose = verbosity
