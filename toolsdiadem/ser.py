import numpy as np

def fetchFreqAndS11( prefix, encoding="ISO-8859-1", verbose=0 ):
    tmp_fre = np.loadtxt( prefix + ".FRE", skiprows=9, encoding=encoding)
    tmp_frp = np.loadtxt( prefix + ".FRP", skiprows=9, encoding=encoding)
    freq = tmp_fre[:,0]
    S11  = ( 10 ** ( 0.05 * tmp_fre[:,1]) ) * np.exp( 1j * np.pi * tmp_frp[:,1] / 180.0 )
    if verbose:
        print(f"freq {freq.shape}, S11 {S11.shape}")
    return freq, S11

def getFreqAndS11( data_dir, seq0=0, verbose=0 ):
    S11_all = np.array([])
    freq_all = np.array([])
    for k in range(3):
        prefix = f"{data_dir}/SEQ_{k+seq0}/S11_ssbande_{k+1}_000"
        if verbose:
            print( prefix )
        seq = k
        ssbande = k + 1
        freq, S11 = fetchFreqAndS11( prefix, verbose=verbose )
        freq_all = np.concatenate((freq_all, freq)) if freq_all.size else freq
        S11_all  = np.concatenate((S11_all,  S11))  if  S11_all.size else S11
    if verbose:
        print(f"freq_all {freq_all.shape}, S11_all {S11_all.shape}")
    return freq_all, S11_all

def getWindow(d, center, width, id="llc"):
    gate = 0 * d + 1e-6
    gate[ np.where( np.abs( d - center ) < width ) ] = 1.0

    d_centre = np.sum( gate * d ) / np.sum( gate )

    gate_a0 = ( 1e-6 + ( ( d - d_centre ) / d_centre ) ** 4 )
    gate_a = 1 / gate_a0

    offset = np.amax( gate_a * ( 1 - gate ) )
    gate_b0 = gate_a / offset - 1.0
    gate_b  = gate_b0 * ( 1 - gate ) + 1.0

    tmp = np.where( gate == 1 )[0]
    y1 = d * 0.0
    y1[ 0 : tmp[0] ] = d[ 0 : tmp[0] ]
    y2 = d * 0.0
    y2[ tmp[-1] : ] = d[ tmp[-1] : ]

    y1 = y1 / d[ tmp[0] ]

    y2[ tmp[-1] : ] = ( y2[ tmp[-1] : ] - np.amax(d) ) / ( d[ tmp[-1] ] - np.amax(d) )

    y1[ tmp[0] ] = 0.0
    y2[ tmp[-1] ] = 0.0

    gate_c = gate_b * ( y1 + y2 + gate )
    gate = gate_c

    return gate