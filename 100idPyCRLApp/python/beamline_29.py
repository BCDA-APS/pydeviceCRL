 
"""
IEX (29ID) Shadow beamline -- forgot which branch this is meant to emulate.

M. Wyman 2021-08-02
"""

__all__ = """
    IEX
""".split()

import Shadow
import numpy as np
import random
import csv

class beamline:

    def __init__(self, elements = 6, dofs = 6):
        self.elementsN = elements
        self.elementsDOFn = dofs

        #pos holds the positions/angles of each element
        self.pos = np.zeros((elements, dofs))
        
        #sRange holds the limits for each element
        self.sRange = np.zeros((elements, dofs, 2))

        #mask is matrix to so that DOFs that aren't used are skipped in the sampling/zeroing/centering process
        self.mask = np.zeros((elements, dofs))
        
        self.resetBeamline()
#        self.beam = Shadow.Beam()
#        self.oe0 = Shadow.Source()
#        self.oe = []
#        for i in range(elements):
#            self.oe.append(Shadow.OE()) 
            
    def zero(self):
        for i in range(self.elementsN):
            for j in range(self.elementsDOFn):
                if (self.mask[i,j] == 1):
                    self.pos[i,j] = 0
#        self.adjust()
        return
        
    def center(self):
        for i in range(self.elementsN):
            for j in range(self.elementsDOFn):
                if (self.mask[i,j] == 1):
                    self.pos[i,j] = float(format((self.sRange[i,j,0]+self.sRange[i,j,1])/2.0))
#        self.adjust()
        return
    
    def sample(self, method = "uniform"):
        for i in range(self.elementsN):
            for j in range(self.elementsDOFn):
                if (self.mask[i,j] == 1):
                    self.pos[i,j] = float(format(random.uniform(self.sRange[i,j,0],self.sRange[i,j,1])))
#        self.adjust()
        return self.pos
          
    
    def adjust(self):
        self.oe[0].OFFX = self.pos[0,0]
#        self.oe[0].OFFY = self.pos[0,1]
#        self.oe[0].OFFZ = self.pos[0,2]    
        self.oe[0].X_ROT = self.pos[0,3]
#        self.oe[0].Y_ROT = self.pos[0,4]
#        self.oe[0].Z_ROT = self.pos[0,5]

        self.oe[1].OFFX = self.pos[1,0]
#        self.oe[1].OFFY = self.pos[1,1]
#        self.oe[1].OFFZ = self.pos[1,2]
        self.oe[1].X_ROT = self.pos[1,3]
#        self.oe[1].Y_ROT = self.pos[1,4]
#        self.oe[1].Z_ROT = self.pos[1,5]

#        self.oe[2].XOFF = self.pos[2,0]
#        self.oe[2].YOFF = self.pos[2,1]
#        self.oe[2].ZOFF = self.pos[2,2]
        self.oe[2].X_ROT = self.pos[2,3]
#        self.oe[2].Y_ROT = self.pos[2,4]
#        self.oe[2].Z_ROT = self.pos[2,5]

#        self.oe[3].XOFF = self.pos[3,0]
#        self.oe[3].YOFF = self.pos[3,1]
#        self.oe[3].ZOFF = self.pos[3,2]
        self.oe[3].X_ROT = self.pos[3,3]
#        self.oe[3].Y_ROT = self.pos[3,4]
#        self.oe[3].Z_ROT = self.pos[3,5]

#        self.oe[4].XOFF = self.pos[4,0]
        self.oe[4].YOFF = self.pos[4,1]
#        self.oe[4].ZOFF = self.pos[4,2]
        self.oe[4].X_ROT = self.pos[4,3]
        self.oe[4].Y_ROT = self.pos[4,4]
#        self.oe[4].Z_ROT = self.pos[4,5]

#        self.oe[5].XOFF = self.pos[5,0]
        self.oe[5].YOFF = self.pos[5,1]
#        self.oe[5].ZOFF = self.pos[5,2]
        self.oe[5].X_ROT = self.pos[5,3]
        self.oe[5].Y_ROT = self.pos[5,4]
#        self.oe[5].Z_ROT = self.pos[5,5]
        return

    def resetBeamline(self):
        self.beam = Shadow.Beam()
        self.oe0 = Shadow.Source()

        self.oe = []
        for i in range(self.elementsN):
            self.oe.append(Shadow.OE())
        
    # Initialize beamline setup
    #  undulator
        self.oe0.FDISTR = 3
        self.oe0.F_COLOR = 3
        self.oe0.F_PHOT = 0
        self.oe0.HDIV1 = 0.5
        self.oe0.HDIV2 = 0.5
        self.oe0.IDO_VX = 0
        self.oe0.IDO_VZ = 0
        self.oe0.IDO_X_S = 0
        self.oe0.IDO_Y_S = 0
        self.oe0.IDO_Z_S = 0
        self.oe0.NPOINT = 100000
        self.oe0.NTOTALPOINT = 20000
        self.oe0.PH1 = 499.97
        self.oe0.PH2 = 500.03
        self.oe0.SIGDIX = 1.964e-05
        self.oe0.SIGDIZ = 1.6349401e-05
        self.oe0.SIGMAX = 0.27609399
        self.oe0.SIGMAZ = 0.0261531994
        self.oe0.VDIV1 = 0.5
        self.oe0.VDIV2 = 0.5

    #  plane mirror
        self.oe[0].ALPHA = 90.0
        self.oe[0].DUMMY = 0.1
        self.oe[0].FHIT_C = 1
        self.oe[0].FWRITE = 3
        self.oe[0].F_MOVE = 1
        self.oe[0].OFFX = self.pos[0,0]
        self.oe[0].OFFY = self.pos[0,1]
        self.oe[0].OFFZ = self.pos[0,2]
        self.oe[0].RLEN1 = 260.0
        self.oe[0].RLEN2 = 260.0
        self.oe[0].RWIDX1 = 20.0
        self.oe[0].RWIDX2 = 20.0
        self.oe[0].T_IMAGE = 0.0
        self.oe[0].T_INCIDENCE = 89.599998
        self.oe[0].T_REFLECTION = 89.599998
        self.oe[0].T_SOURCE = 30800.0
        self.oe[0].X_ROT = self.pos[0,3]
        self.oe[0].Y_ROT = self.pos[0,4]
        self.oe[0].Z_ROT = self.pos[0,5]
    #  when taking into account reflectivity
    # oe1.FILE_REFL = b"C:/cygwin64/Oasys/Si.dat"
    # oe1.F_REFLEC = 1

    #  plane mirror
        self.oe[1].DUMMY = 0.1
        self.oe[1].FHIT_C = 1
        self.oe[1].FWRITE = 3
        self.oe[1].F_MOVE = 1
        self.oe[1].OFFX = self.pos[1,0]
        self.oe[1].OFFY = self.pos[1,1]
        self.oe[1].OFFZ = self.pos[1,2]
        self.oe[1].RLEN1 = 75.0
        self.oe[1].RLEN2 = 75.0
        self.oe[1].RWIDX1 = 22.5
        self.oe[1].RWIDX2 = 22.5
        self.oe[1].T_IMAGE = 0.0
        self.oe[1].T_INCIDENCE = 88.5
        self.oe[1].T_REFLECTION = 88.5
        self.oe[1].T_SOURCE = 500.0
        self.oe[1].X_ROT = self.pos[1,3]
        self.oe[1].Y_ROT = self.pos[1,4]
        self.oe[1].Z_ROT = self.pos[1,5]
    #  reflect
    # oe2.FILE_REFL = b"C:/cygwin64/Oasys/Si.dat"
    # oe2.F_REFLEC = 1

    #  plane mirror
        self.oe[2].ALPHA = 90.0
        self.oe[2].DUMMY = 0.1
        self.oe[2].FHIT_C = 1
        self.oe[2].FWRITE = 3
        self.oe[2].F_MOVE = 1
        self.oe[2].XOFF = self.pos[2,0]
        self.oe[2].YOFF = self.pos[2,1]
        self.oe[2].ZOFF = self.pos[2,2]
        self.oe[2].RLEN1 = 190.0
        self.oe[2].RLEN2 = 190.0
        self.oe[2].RWIDX1 = 15.0
        self.oe[2].RWIDX2 = 15.0
        self.oe[2].T_IMAGE = 0.0
        self.oe[2].T_INCIDENCE = 86.810772
        self.oe[2].T_REFLECTION = 86.810772
        self.oe[2].T_SOURCE = 8265.8165
        self.oe[2].X_ROT = self.pos[2,3]
        self.oe[2].Y_ROT = self.pos[2,4]
        self.oe[2].Z_ROT = self.pos[2,5]
    #  reflect
    # oe3.FILE_REFL = b"C:/cygwin64/Oasys/Si.dat"
    # oe3.F_REFLEC = 1

    #  plane grating
        self.oe[3].ALPHA = 180.0
        self.oe[3].DUMMY = 0.1
        self.oe[3].FHIT_C = 1
        self.oe[3].FWRITE = 2
        self.oe[3].F_GRATING = 1
        self.oe[3].F_MOVE = 1
        self.oe[3].F_RULING = 5
        self.oe[3].F_RUL_ABS = 1
        self.oe[3].XOFF = self.pos[3,0]
        self.oe[3].YOFF = self.pos[3,1]
        self.oe[3].ZOFF = self.pos[3,2]
        self.oe[3].RLEN1 = 57.5
        self.oe[3].RLEN2 = 57.5
        self.oe[3].RULING = 1199.22002
        self.oe[3].RUL_A1 = 0.165491998
        self.oe[3].RUL_A2 = 1.0793e-05
        self.oe[3].RUL_A3 = 1.99999999e-06
        self.oe[3].RWIDX1 = 12.5
        self.oe[3].RWIDX2 = 12.5
        self.oe[3].T_IMAGE = 2000.0
        self.oe[3].T_INCIDENCE = 87.771697
        self.oe[3].T_REFLECTION = 85.049847
        self.oe[3].T_SOURCE = 135.0
        self.oe[3].X_ROT = self.pos[3,3]
        self.oe[3].Y_ROT = self.pos[3,4]
        self.oe[3].Z_ROT = self.pos[3,5]

    # ellipsoid mirror oe6 -> oe5 after removing slit
        self.oe[4].ALPHA = 90.0
        self.oe[4].AXMAJ = 33150.0
        self.oe[4].AXMIN = 296.85199
        self.oe[4].DUMMY = 0.1
        self.oe[4].ELL_THE = 0.186719999
        self.oe[4].FCYL = 1
        self.oe[4].FHIT_C = 1
        self.oe[4].FMIRR = 2
        self.oe[4].FWRITE = 3
        self.oe[4].F_EXT = 1
        self.oe[4].F_MOVE = 1
        self.oe[4].XOFF = self.pos[4,0]
        self.oe[4].YOFF = self.pos[4,1]
        self.oe[4].ZOFF = self.pos[4,2]
        self.oe[4].RLEN1 = 140.0
        self.oe[4].RLEN2 = 140.0
        self.oe[4].RWIDX1 = 7.5
        self.oe[4].RWIDX2 = 7.5
        self.oe[4].T_IMAGE = 0.0
        self.oe[4].T_INCIDENCE = 88.5
        self.oe[4].T_REFLECTION = 88.5
        self.oe[4].T_SOURCE = 4600.0
        self.oe[4].X_ROT = self.pos[4,3]
        self.oe[4].Y_ROT = self.pos[4,4]
        self.oe[4].Z_ROT = self.pos[4,5]
    #  reflect
    # oe6.FILE_REFL = b"C:/cygwin64/Oasys/Si.dat"
    # oe6.F_REFLEC = 1

    #  ellipsoid mirror oe7 -> oe6 after removing slit
        self.oe[5].ALPHA = 90.0
        self.oe[5].AXMAJ = 3300.0
        self.oe[5].AXMIN = 66.6355972
        self.oe[5].DUMMY = 0.1
        self.oe[5].ELL_THE = 1.40167999
        self.oe[5].FCYL = 1
        self.oe[5].FHIT_C = 1
        self.oe[5].FMIRR = 2
        self.oe[5].FWRITE = 3
        self.oe[5].F_EXT = 1
        self.oe[5].F_MOVE = 1
        self.oe[5].XOFF = self.pos[5,0]
        self.oe[5].YOFF = self.pos[5,1]
        self.oe[5].ZOFF = self.pos[5,2]
        self.oe[5].RLEN1 = 35.0
        self.oe[5].RLEN2 = 35.0
        self.oe[5].RWIDX1 = 10.0
        self.oe[5].RWIDX2 = 10.0
        self.oe[5].T_IMAGE = 1200.0
        self.oe[5].T_INCIDENCE = 88.5
        self.oe[5].T_REFLECTION = 88.5
        self.oe[5].T_SOURCE = 800.0
        self.oe[5].X_ROT = self.pos[5,3]
        self.oe[5].Y_ROT = self.pos[5,4]
        self.oe[5].Z_ROT = self.pos[5,5]  
        
        return
    
    def run(self, iwrite = False, verbose = False, **kwargs):
        self.resetBeamline()        
        self.beam.genSource(self.oe0)

        if iwrite:
            self.oe0.write("end.00")
            self.beam.write("begin.dat")
            
        for i in range(6):
            if verbose: print("\r    Running optical element: %d" % (i+1), end = "")

            if iwrite:
                self.oe[i].write("start"+'0'+str(i+1))
            self.beam.traceOE(self.oe[i], i+1)
            if iwrite:
                self.oe[i].write("end."+'0'+str(i+1))
                self.beam.write("star."+'0'+str(i+1))

        if kwargs:
            results = self.beam.histo2(1, 2, **kwargs)
        else:
            results = self.beam.nrays(nolost=1)
        return results
    
IEX = beamline(elements = 6)




IEX.sRange = np.asarray([[[-20.4, 20.4],[-20., 20.],[-1.5, 1.5],[-0.01, 0.01],[-0.62, 0.62],[-20, 20]],
                         [[-22.9, 22.9],[-20., 20.],[-2.5, 2.5],[-0.01, 0.01],[-0.18, 0.18],[-20, 20]],
                         [[-16.1, 16.1],[-20., 20.],[-0.9, 0.9],[-0.02, 0.02],[-0.37, 0.37],[-20, 20]],
                         [[-13.6, 13.6],[-3.5, 3.5],[-1., 1.],[-0.02, 0.02],[-0.35, 0.35],[-7.49, 7.49]],
                         [[-8.3, 8.3],[-42, 42],[-5., 5.],[-0.4, 0.4],[-3.22, 3.22],[-20, 20]],
                         [[-10.8, 10.8],[-86.3, 86.3],[-2.3, 2.3],[-1.38, 1.38],[-8.7, 8.7],[-20, 20]]])

#IEX (forgot which branch) really only has 12 degrees of freedom
IEX.mask = np.asarray([[1, 0, 0, 1, 0, 0],
                       [1, 0, 0, 1, 0, 0],
                       [0, 0, 0, 1, 0, 0],
                       [0, 0, 0, 1, 0, 0], 
                       [0, 1, 0, 1, 1, 0],
                       [0, 1, 0, 1, 1, 0]])
