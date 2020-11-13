import os

import numpy as np

def readCfg( cfg_file ):
    
    CFG = {}
    pi = np.pi
    
    with open( cfg_file ) as file:
        for k in range(4):
            file.readline()
        H_nom = file.readline().split(':')[1]
        H_nom = H_nom.replace( ',', '.' )
        CFG["H_nom"]     = float(H_nom.strip())
        CFG["xamin_nom"] = float(file.readline().split(':')[1])
        CFG["xamax_nom"] = float(file.readline().split(':')[1])
        CFG["Npos"]      = int(float(file.readline().rsplit(':')[1]))
        CFG["theta_c"]   = float(file.readline().split(':')[1]) * pi / 180
        CFG["theta_ap"]  = float(file.readline().split(':')[1]) * pi / 180
        CFG["phi_c"]     = float(file.readline().split(':')[1]) * pi / 180
        CFG["phi_ap"]    = float(file.readline().split(':')[1]) * pi / 180
        for k in range(3):
            file.readline()
        CFG["Fmin"] = float(file.readline().split(':')[1])
        CFG["Fmax"] = float(file.readline().split(':')[1])
        CFG["Nf"]   = int(float(file.readline().split(':')[1]))
        
        CFG["c"] = 3e8
    
    return CFG

def readFile( filename, Npos, Nf ):
    # pos0f0 .. pos0fNf pos1f0 .. pos1fNf posNposf0 .. posNposfNf
    with open( filename ) as file:
        dum = np.fromfile(file, dtype = np.float32)
        c = (dum[ 0 : : 2 ] + 1j * dum[ 1 : : 2 ]).reshape(Npos, Nf).astype(complex)
    return c

############
############
## POSITIONS

class AntPosPoSARGB(object):
    """docstring for AntPosT"""
    def __init__(self, Txx_Txy_Txz_Rxx_Rxy_Rxz):
        super(AntPosPoSARGB, self).__init__()

        self.Txx = Txx_Txy_Txz_Rxx_Rxy_Rxz[0,:]
        self.Txy = Txx_Txy_Txz_Rxx_Rxy_Rxz[1,:]
        self.Txz = Txx_Txy_Txz_Rxx_Rxy_Rxz[2,:]
        self.Rxx = Txx_Txy_Txz_Rxx_Rxy_Rxz[3,:]
        self.Rxy = Txx_Txy_Txz_Rxx_Rxy_Rxz[4,:]
        self.Rxz = Txx_Txy_Txz_Rxx_Rxy_Rxz[5,:]
        self.mean_x = ( self.Txx + self.Rxx ) / 2
        self.mean_y = ( self.Txy + self.Rxy ) / 2
        self.mean_z = ( self.Txz + self.Rxz ) / 2

def readPos( filename ):
    with open( filename ) as file:
        ant_pos = AntPosPoSARGB( np.fromfile(file, dtype = np.float32).reshape((-1, 6)).T )
        return ant_pos

def getZRail( s, units="" ):
    if units=="cm":
        zRail = float( s.split('cm')[0] ) / 100
    else:
        zRail = float( s.split('m')[0] + '.' + s.split('m')[1] )
    return zRail

def buildPos(xMin, xMax, Nx, Dx1, Dx2, h, zRef):
    offsetX = (Dx1+2*Dx2) / 2
    offsetY = 0
    offsetZ = 0
    A1 = np.array([ 0,           0, 4 * h ] )
    A2 = np.array([ Dx1,         0, 2 * h ] )
    A3 = np.array([ Dx1+Dx2,     0, 1 * h ] )
    A4 = np.array([ Dx1+Dx2+Dx2, 0, 0 * h ] )
    col0 = np.zeros( Nx ).T.reshape(-1,1)
    col1 = np.arange( Nx ).T.reshape(-1,1)
    xTrack = np.linspace( xMin, xMax, Nx ) - offsetX
    yTrack = np.zeros( Nx )                - offsetY
    zTrack = np.ones( Nx ) * zRef          - offsetZ
    track = np.array([xTrack, yTrack, zTrack]).T
    xyzAnt = {}
    xyzAnt["1"] = np.hstack( ( col0, col1, A1 + track ) )
    xyzAnt["2"] = np.hstack( ( col0, col1, A2 + track ) )
    xyzAnt["3"] = np.hstack( ( col0, col1, A3 + track ) )
    xyzAnt["4"] = np.hstack( ( col0, col1, A4 + track ) )
    return xyzAnt

def buildPosUsingTrack(prefix, Dx1, Dx2, h, zRef, verbose=0, old=0, bizona=(0,0,0)):
    A1 = np.array([ 0,           0, 4 * h ] )
    A2 = np.array([ Dx1,         0, 2 * h ] )
    A3 = np.array([ Dx1+Dx2,     0, 1 * h ] )
    A4 = np.array([ Dx1+Dx2+Dx2, 0, 0 * h ] )
    if old:
        filename = prefix + "21.dat.old"
    else:
        filename = prefix + "21.dat"
    # antenna 1
    with open( filename ) as file:
        track = np.fromfile(file, dtype = np.float32).reshape(-1, 6)
        xTrack = track[:,0].reshape(-1,1)
        Npos = xTrack.size
        yTrack = np.zeros( xTrack.shape )
        zTrack = np.ones(  xTrack.shape ) * zRef
        track = np.hstack([xTrack, yTrack, zTrack])
        col0 = np.zeros( xTrack.shape )
        col1 = np.arange( Npos ).reshape(-1,1)
        xyzAnt = {}

        if bizona == (0,0,0,):
            xyzAnt["1"] = np.hstack( ( col0, col1, A1 + track ) )
            xyzAnt["2"] = np.hstack( ( col0, col1, A2 + track ) )
            xyzAnt["3"] = np.hstack( ( col0, col1, A3 + track ) )
            xyzAnt["4"] = np.hstack( ( col0, col1, A4 + track ) )

        if bizona != (0,0,0):
            xyzAnt["1"] = np.hstack( ( col0, col1, A2 + track ) )
            xyzAnt["2"] = np.hstack( ( col0, col1, A3 + track ) )
            xyzAnt["3"] = np.hstack( ( col0, col1, A4 + track ) )
            # 4
            xyzAnt["4"] = np.hstack( ( col0, col1, A4 + track ) )
            xyzAnt["4"][:,2] = np.ones( Npos ) * bizona[0]
            xyzAnt["4"][:,3] = np.ones( Npos ) * bizona[1]
            xyzAnt["4"][:,4] = np.ones( Npos ) * bizona[2]

    if verbose:
        print( f"read {Npos} positions in file:\n{filename}" )
    
    return xyzAnt
   
######
######
## CAL

def readCal( cal_file ):
    
    CAL = {}
    
    with open( cal_file ) as file:
        for k in range(4):
            file.readline()
        
        CAL["xc"]   = float(file.readline().split(':')[1])
        CAL["rc"]   = float(file.readline().split(':')[1])
        CAL["rmse"] = float(file.readline().split(':')[1])
    
        for k in range(3):
            file.readline()
        
        CAL["ref_rg_bin"] = float(file.readline().split(':')[1])
        CAL["d_shift"]    = float(file.readline().split(':')[1])

        return CAL

def readData( data_path, img_name, cal_name="", verbose=0 ):
    
    # Read configuration
    CFG = readCfg( data_path + "/PoSAR.cfg" )
    
    # Read raw data
    RD = readFile( data_path + "/" + img_name, CFG['Npos'], CFG['Nf'] )
    maxRD = np.amax( abs(RD) )
    if verbose:
        print( f"maxRD = {maxRD}" )
    
    # Read calibration information
    CAL = {}
    if cal_name:
        print("WARNING cal_name parameter is not empty, one should check the behaviour!")
        cal_path = data_path
        img_name_without_ext = os.path.splitext( img_name )[0]
        CAL = readCal( cal_name )
    else:
        CAL["d_shift"] = 0 # 0.2
    
    # Read frequency coordinates
    with open( data_path + "/freq.dat" ) as file:
        f = np.fromfile( file, dtype = np.float32 )
    
    Fmin = min( f )
    Fmax = max( f )
    B = Fmax - Fmin
    delta_f = B / (CFG["Nf"]-1)
    f_kc = f[ int((CFG["Nf"]+1)/2 ) ]
    kc = 4 * np.pi * f_kc / CFG["c"]
    CFG["Fmin"] = Fmin
    CFG["Fmax"] = Fmax
    CFG["kc"] = kc
    
    if verbose:
        print(f"Fmin = {Fmin:.3e}, Fmax = {Fmax:.3e}, B = {B:.3e}")
        print(f"kc = {kc:.3f} for f = {f_kc:.3e}, delta_f = {delta_f:e}")
    
    CFG_Fmin = CFG["Fmin"]
    CFG_Fmax = CFG["Fmax"]
    CFG_B = CFG_Fmax - CFG_Fmin
    if verbose:
        print(f"CFG_Fmin = {CFG_Fmin:.3e}, Fmax = {CFG_Fmax:.3e}, CFG_B = {CFG_B:.3e}")
    
    delta_sr = CFG["c"] / (2 * B)
    CFG["delta_sr"] = delta_sr
    if verbose:
        print(f"delta_sr = {delta_sr:.3e}")
    
    return CFG, RD, CAL, f

