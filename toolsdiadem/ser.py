import numpy as np
import re, os

def correctBin(S11, fBin):
    S11[:,fBin] = (S11[:,fBin+1] + S11[:,fBin-1]) / 2

def computeF2C(base_path, dic_dir, nb_elev, nb_freq, nb_ssb, tag=""):
    nbAbs = 0
    for key in dic_dir:
        nbAbs = nbAbs + len(dic_dir[key])
    print( f"nbAbs = {nbAbs}" )
    S11_moy = np.array([])
    for key in dic_dir:
        print(key)
        for absx in dic_dir[key]:
            print(absx)
            data_dir = f"{base_path}/{key}/{absx}" 
            freq_S11, S11 = getData_ssb(nb_elev, nb_freq, nb_ssb, data_dir, tag=tag)
            if S11_moy.shape == (0,):
                S11_moy = np.mean(S11, axis=0)
            else:
                S11_moy += np.mean(S11, axis=0)
            
    S11_moy /= nbAbs
    return S11_moy

def dB(x):
    y = 20 * np.log10( np.abs( x ) )
    return y

def getData_ssb( nb_elev, nb_freq, nb_ssb, data_dir, encoding="ISO-8859-1", tag=""):
    incr = 0
    S11 = np.zeros( (nb_elev, nb_freq), dtype=complex )
    freq_elev = np.zeros( nb_freq )
    S11_elev = np.zeros( nb_freq, dtype=complex )

    for ielev in range( nb_elev ):
        counter = 0
        for ssb in range( nb_ssb ):
            prefix = f"{data_dir}/SEQ_{incr}/{tag}S11_ssbande_{ssb+1}_000"
            tmp_fre = np.loadtxt( prefix + ".FRE", skiprows=9, encoding=encoding)
            tmp_frp = np.loadtxt( prefix + ".FRP", skiprows=9, encoding=encoding)
            freq_ssb = tmp_fre[:,0]
            S11_ssb  = ( 10 ** ( 0.05 * tmp_fre[:,1]) ) * np.exp( 1j * np.pi * tmp_frp[:,1] / 180.0 )
            nb_freq_ssb = freq_ssb.size
            freq_elev[ counter : counter + nb_freq_ssb ] = freq_ssb
            S11_elev[  counter : counter + nb_freq_ssb ] = S11_ssb
            counter = counter + nb_freq_ssb
            incr = incr + 1

        S11[ ielev, : ] = S11_elev

    return freq_elev, S11

def getAzSlRoElPolFFPolCATR( seqPath, encoding="ISO-8859-1" ):
    nbSEQ = 0
    p = re.compile("^[0-9]+$")
    for root, dirs, files in os.walk(seqPath):
        number = root.split("/")[-1].split("SEQ_")[-1]
        if p.match(number) == None:
            nbSEQ += 1
    azSlRoElPolFFPolCATR = np.zeros((nbSEQ, 6))
    for seq in range(nbSEQ):
        filename = seqPath + f"/SEQ_{seq}/Axis.pos"
        vals = np.loadtxt( filename, delimiter="=", usecols=1)
        azSlRoElPolFFPolCATR[seq, :] = vals
    return azSlRoElPolFFPolCATR

def fetchFreqAndS11( prefix, encoding="ISO-8859-1", verbose=0 ):
    tmp_fre = np.loadtxt( prefix + ".FRE", skiprows=9, encoding=encoding)
    tmp_frp = np.loadtxt( prefix + ".FRP", skiprows=9, encoding=encoding)
    freq = tmp_fre[:,0]
    S11  = ( 10 ** ( 0.05 * tmp_fre[:,1]) ) * np.exp( 1j * np.pi * tmp_frp[:,1] / 180.0 )
    if verbose:
        print(f"freq {freq.shape}, S11 {S11.shape}")
    return freq, S11

def getFreqAndS11( data_dir, nb_ssb, seq0=0, verbose=0 ):
    S11_all = np.array([])
    freq_all = np.array([])
    for k in range(nb_ssb):
        prefix = f"{data_dir}/SEQ_{k+seq0}/S11_ssbande_{k+1}_000"
        if verbose:
            print( prefix )
        freq, S11 = fetchFreqAndS11(prefix, verbose=verbose)
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

def getReflectivite( S11, plaqueRef, fondDeChambre, gate ):
    S11_plaque_ref_tg = np.fft.fft(np.fft.ifft(plaqueRef - fondDeChambre) * gate)
    S11_minus_fdc    = S11 - fondDeChambre
    s11_minus_fdc_td = np.fft.ifft(S11_minus_fdc, axis=1 )
    s11_minus_fdc_tg = s11_minus_fdc_td * gate
    S11_minus_fdc_tg = np.fft.fft(s11_minus_fdc_tg, axis=1 )
    reflectivite     = S11_minus_fdc_tg / S11_plaque_ref_tg
    return reflectivite