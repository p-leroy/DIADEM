import os
import fnmatch

import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1 import AxesGrid
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt

def get_Nf_Nstep(data_dir):
    f_step = []
    for file in os.listdir(data_dir):
        if fnmatch.fnmatch(file, "*.amp"):
            sfile = file.split(".amp")[0].split("_")
            f_step.append((sfile[1], sfile[2]))
    Nf = max([int(f_step[0]) for f_step in f_step if f_step[1] == '000']) + 1
    Nstep = max([int(f_step[1]) for f_step in f_step if f_step[0] == '000']) + 1
    print(f"data_dir: {data_dir}")
    print("Nf = {}, Nstep = {}".format(Nf, Nstep))
    return Nf, Nstep

def getF(nbFreq, data_prefix, root_dir):
    f = np.zeros(nbFreq)
    for freqIdx in range(nbFreq):
        filename = os.path.join(root_dir, f"{data_prefix}_{freqIdx:03d}_000.amp")
        with open(filename) as file:
            for line in file:
                if "FreqValue" in line:
                    f[freqIdx] = float(line.split("=")[1])
                elif "Data#1" in line:
                    break
    print(f"data_dir: {root_dir}")
    print("Nf = {}, fmin = {}, fmax = {}".format(f.size, np.amin(f), np.amax(f)))
    return f

def getStep(nbY, data_prefix, root_dir):
    step = np.zeros(nbY)
    for idx in range(nbY):
        filename = os.path.join(root_dir, f"{data_prefix}_000_{idx:03d}.amp")
        with open(filename) as file:
            for line in file:
                if "StepAxis=" in line:
                    step[idx] = float(line.split(" ")[-2])
                    axis = line.split(" ")[0].split("=")[1]
                elif "Data#1" in line:
                    break
    print(f"data_dir: {root_dir}")
    print(f"Nstep {step.size} ({axis}), min {np.amin(step)}, max {np.amax(step)}")
    return step

def get_headerSize(filename):
    with open(filename) as file:
        for counter, line in enumerate(file):
            if "Data#1" in line:
                break
    return counter + 1

def get_axes(filename):
    with open(filename) as file:
        for line in file:
            if "ScanAxis=" in line:
                ScanAxis = line.split("=")[1].rstrip()
            if "StepAxis=" in line:
                StepAxis = line.split(" ")[0].split("=")[1]
                StepPosition = float(line.split(" ")[-2])
            elif "Data#1" in line:
                break
    return ScanAxis, StepAxis, StepPosition

def get_scan_val(data_prefix, fi, yi, ext):
    filename = data_prefix + f"_{fi:03d}_{yi:03d}" + ext
    headerSize = get_headerSize(filename)
    ScanAxis, StepAxis, StepPosition = get_axes(filename)
    Freq = None
    with open(filename) as file:
        for line in file:
            if "Remarks=Freq." in line:
                freq = float(line.split(" ")[-2])
            elif "Data#1" in line:
                break
    # read the data
    scan, val = np.genfromtxt(filename, skip_header=headerSize, unpack=True)
    print(f"[{ext}] {freq} GHz, Nscan {scan.size} ({ScanAxis}), min {np.amin(scan)}, max {np.amax(scan)}, StepAxis {StepAxis}, StepPosition {StepPosition}")
    return scan, val, freq

def getAmpPhaArrays(Nscan, Nstep, fi, data_prefix):
    ampArray = np.zeros((Nstep, Nscan))
    phaArray = np.zeros((Nstep, Nscan))
    for step in range(Nstep):
        scan, amp, freq = get_scan_val(data_prefix, fi, step, ".amp")
        ampArray[step, :] = amp[:]
        scan, pha, freq = get_scan_val(data_prefix, fi, step, ".pha")
        phaArray[step, :] = pha[:]
    return ampArray, phaArray, freq

def getComplex(amp, pha):
    linAmp = np.power( 10, amp / 20 )
    re = linAmp * np.cos( pha * np.pi / 180 )
    im = linAmp * np.sin( pha * np.pi / 180 )
    return re + 1j * im

def getThetaPhi(kx, ky, kz, k0):
    theta = np.arccos( kz / k0 )
    phi = np.empty(kx.shape)
    phi = np.arctan2(ky, kx)
    return theta, phi

def saveData(root_dir, ampAndPha, scan, freq, step):
    header = "scan amp pha"
    filename = root_dir + f"{freq}_step{step}_amp_pha_32GHz.data"
    print("save " + filename)
    np.savetxt(filename, 
               np.column_stack(( 
                   scan.flatten(),
                   ampAndPha[0][step].flatten(), 
                   ampAndPha[1][step].flatten())),
               header=header,
               comments=""
              )

def addColorBar(im, ax, aspect='equal'):
    ax.grid()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad="5%")
    plt.colorbar(im, cax=cax)
    ax.set_aspect(aspect)

def E_E_yaghjian(theta, k, a, b, AE):
    beta_over_k = (1 - ( np.pi / (k*a) )**2)**0.5 # normalized propagation constant for the TE10 mode
    E_E_theta = np.zeros(theta.shape)
    idx = np.where(theta != 0)
    E_E_theta[idx] = AE * ( 1 + beta_over_k*np.cos(theta[idx]) ) / ( 1 + beta_over_k ) \
                   * np.sin( k*b/2 * np.sin(theta[idx]) ) / ( k*b/2 * np.sin(theta[idx]) )
    idx = np.where(theta == 0)
    E_E_theta[idx] = AE
    return E_E_theta

def E_H_yaghjian(theta, k, a, b, AE):
    E_H_theta = np.zeros(theta.shape)
    E_H_theta = AE * (np.pi/2)**2 * np.cos(theta) \
                   * np.cos(k*a/2*np.sin(theta)) / ((np.pi/2)**2 - (k*a/2*np.sin(theta))**2)
    return E_H_theta

def E_theta_phi_yaghjian(theta, phi, k, a, b, AE, inv=1):
    beta_over_k = (1 - ( np.pi / (k*a) )**2)**0.5 # normalized propagation constant for the TE10 mode

    E_E_theta = AE * ( 1 + beta_over_k*np.cos(theta) ) / ( 1 + beta_over_k ) \
                    * np.sinc( k*b/2 * np.sin(theta) / np.pi )                   

    E_H_theta = AE * (np.pi/2)**2 * np.cos(theta) \
                   * np.cos(k*a/2*np.sin(theta)) / ((np.pi/2)**2 - (k*a/2*np.sin(theta))**2)

    return E_E_theta * np.sin(phi), E_H_theta * np.cos(phi)

def plot_E_LHCP_E_RHCP(u, v, E_LHCP, E_RHCP, nameX, nameY, vmax=-1, vmin=-1, suffix=""):
    fig, (ax0, ax1) = plt.subplots(1, 2)

    if vmin==-1 and vmax ==-1:
        im0 = ax0.pcolor(u, v, 20*np.log10(np.abs(E_LHCP)))
        im1 = ax1.pcolor(u, v, 20*np.log10(np.abs(E_RHCP)))
    else:
        im0 = ax0.pcolor(u, v, 20*np.log10(np.abs(E_LHCP)), vmax=vmax, vmin=vmin)
        im1 = ax1.pcolor(u, v, 20*np.log10(np.abs(E_RHCP)), vmax=vmax, vmin=vmin)


    ax0.set_title("E_LHCP"+suffix)
    ax0.set_xlabel("u")
    ax0.set_ylabel("v")
    addColorBar(im0, ax0)

    ax1.set_title("E_RHCP"+suffix)
    ax1.set_xlabel("u")
    ax1.set_ylabel("v")
    addColorBar(im1, ax1)

    title = "data X: " + nameX + "\ndata Y: " + nameY
    fig.suptitle(title)

def plotImgMosaic( imgOrig, extent, origin='upper', T=0 ):
    fig = plt.figure()
    grid = AxesGrid(fig, 111,  # similar to subplot(111)
                    nrows_ncols=(3, 2),
                    axes_pad=0.0,
                    share_all=True,
                    label_mode="all",
                    cbar_location="top",
                    cbar_mode="single",
                   )
    
    img = dict(imgOrig)

    for im in img:
        img[im] = np.squeeze(img[im])
        if T:
            img[im] = img[im].T
    
    ij = "21"
    toPlot = 20 * np.log10( np.abs( img[ij] ) )
    im = grid[0].imshow(toPlot, cmap='jet', extent=extent, origin=origin)
    grid[0].set_ylabel(ij)

    ij = "41"
    toPlot = 20 * np.log10( np.abs( img[ij] ) )
    im = grid[1].imshow(toPlot, cmap='jet', extent=extent, origin=origin)
    grid[1].set_ylabel(ij)
    grid[1].yaxis.set_ticks_position('right')
    grid[1].yaxis.set_label_position('right')

    ij = "31"
    toPlot = 20 * np.log10( np.abs( img[ij] ) )
    im = grid[2].imshow(toPlot, cmap='jet', extent=extent, origin=origin)
    grid[2].set_ylabel(ij)

    ij = "42"
    toPlot = 20 * np.log10( np.abs( img[ij] ) )
    im = grid[3].imshow(toPlot, cmap='jet', extent=extent, origin=origin)
    grid[3].set_ylabel(ij)
    grid[3].yaxis.set_ticks_position('right')
    grid[3].yaxis.set_label_position('right')

    ij = "32"
    toPlot = 20 * np.log10( np.abs( img[ij] ) )
    im = grid[4].imshow(toPlot, cmap='jet', extent=extent, origin=origin)
    grid[4].set_ylabel(ij)

    ij = "43"
    toPlot = 20 * np.log10( np.abs( img[ij] ) )
    im = grid[5].imshow(toPlot, cmap='jet', extent=extent, origin=origin)
    grid[5].set_ylabel(ij)
    grid[5].yaxis.set_ticks_position('right')
    grid[5].yaxis.set_label_position('right')

    grid.cbar_axes[0].colorbar(im)
    
    return fig, grid

def plotImgMosaicAutoCross( I1, I2, I3, I12, I13, I23, extent, origin='upper', vmin=-1, vmax=-1 ):
    pi = np.pi
    fig = plt.figure()
    grid = AxesGrid(fig, 111,  # similar to subplot(111)
                    nrows_ncols=(3, 3),
                    axes_pad=0.0,
                    share_all=True,
                    label_mode="L",
                    cbar_location="top",
                    cbar_mode="single", # single each edge
                   )

    tmpMin = 10 * np.log10( np.amin( ( I1, I2, I3 ) ) )
    tmpMax = 10 * np.log10( np.amax( ( I1, I2, I3 ) ) )
    print( f"min = {tmpMin:.2f}, max = {tmpMax:.2f}" )
    if vmin==-1:
        vmin = tmpMin
        vmax = tmpMax
    
    ij = "11"
    toPlot = 10 * np.log10( I1 )
    imI = grid[0].imshow(toPlot, origin='lower', cmap='jet', extent=extent, vmin=vmin,vmax=vmax)
    grid[0].set_ylabel(ij)

    ij = "12"
    toPlot = np.abs( I12 )
    imCross = grid[1].imshow(toPlot, origin='lower', cmap='gray', extent=extent, vmin=0, vmax=1)

    ij = "13"
    toPlot = np.abs( I13 )
    im = grid[2].imshow(toPlot, origin='lower', cmap='gray', extent=extent, vmin=0, vmax=1)
    grid[2].set_ylabel(ij)
    grid[2].yaxis.set_ticks_position('right')
    grid[2].yaxis.set_label_position('right')

    ij = "12"
    toPlot = np.angle( I12 )
    imAngle = grid[3].imshow(toPlot, origin='lower', cmap='jet', extent=extent, vmin=-pi, vmax=pi)
    grid[3].set_ylabel(ij)

    ij = "22"
    toPlot = 10 * np.log10( I1 )
    im = grid[4].imshow(toPlot, origin='lower', cmap='jet', extent=extent, vmin=vmin,vmax=vmax)

    ij = "23"
    toPlot = np.abs( I23 )
    im = grid[5].imshow(toPlot, origin='lower', cmap='gray', extent=extent, vmin=0, vmax=1)
    grid[5].set_ylabel(ij)
    grid[5].yaxis.set_ticks_position('right')
    grid[5].yaxis.set_label_position('right')

    ij = "13"
    toPlot = np.angle( I13 )
    im = grid[6].imshow(toPlot, origin='lower', cmap='jet', extent=extent, vmin=-pi, vmax=pi)
    grid[6].set_ylabel(ij)
    
    ij = "23"
    toPlot = np.angle( I23 )
    im = grid[7].imshow(toPlot, origin='lower', cmap='jet', extent=extent, vmin=-pi, vmax=pi)

    ij = "33"
    toPlot = 10 * np.log10( I3 )
    im = grid[8].imshow(toPlot, origin='lower', cmap='jet', extent=extent, vmin=vmin,vmax=vmax)
    grid[8].set_ylabel(ij)
    grid[8].yaxis.set_ticks_position('right')
    grid[8].yaxis.set_label_position('right')

    grid.cbar_axes[0].colorbar(imI)

    for ax in grid:
        ax.xaxis.set_major_locator(MaxNLocator(prune='both', integer=True))
        ax.yaxis.set_major_locator(MaxNLocator(prune='both', integer=True))
    
    return fig, grid

def plotImgMosaicMonoBi( imgMono, imgBi, extent, origin='upper', T=0, vmin=-1,vmax=-1 ):
    fig = plt.figure()
    grid = AxesGrid(fig, 111,  # similar to subplot(111)
                    nrows_ncols=(2, 1),
                    axes_pad=0.0,
                    share_all=True,
                    label_mode="all",
                    cbar_location="top",
                    cbar_mode="single",
                   )

    if vmin == -1:
        toPlot = 20 * np.log10( np.abs( imgMono ) )
        im = grid[0].imshow(toPlot, cmap='jet', extent=extent, origin=origin)
        grid[0].set_ylabel("monostatic")
        toPlot = 20 * np.log10( np.abs( imgBi ) )
        im = grid[1].imshow(toPlot, cmap='jet', extent=extent, origin=origin)
        grid[1].set_ylabel("bistatic")
    else:
        toPlot = 20 * np.log10( np.abs( imgMono ) )
        im = grid[0].imshow(toPlot, cmap='jet', extent=extent, origin=origin, vmin=vmin, vmax=vmax)
        grid[0].set_ylabel("monostatic")
        toPlot = 20 * np.log10( np.abs( imgBi ) )
        im = grid[1].imshow(toPlot, cmap='jet', extent=extent, origin=origin, vmin=vmin, vmax=vmax)
        grid[1].set_ylabel("bistatic")

    grid.cbar_axes[0].colorbar(im)
    
    return fig, grid

def plotImgMosaicMonoBi2DTomo( mono2D, bi2D, monoTomo, biTomo, extentyx, extentyz, vmin=-1,vmax=-1 ):
    fig = plt.figure()
    grid0 = ImageGrid(fig, 211,  # similar to subplot(111)
                    nrows_ncols=(1, 2),
                    axes_pad=0.0,
                    share_all=True,
                    label_mode="all",
                    cbar_location="top",
                    cbar_mode="single",
                   )
    grid1 = ImageGrid(fig, 212,  # similar to subplot(111)
                    nrows_ncols=(1, 2),
                    axes_pad=0.0,
                    share_all=True,
                    label_mode="all",
                    cbar_location="top",
                    cbar_mode="single",
                   )

    max2D   = np.amax( np.abs( ( mono2D, bi2D ) ) )
    maxTomo = np.amax( np.abs( ( monoTomo, biTomo ) ) )
    vmaxTmp    = 20 * np.log10( max(max2D, maxTomo) )
    min2D   = np.amin( np.abs( ( mono2D, bi2D ) ) )
    minTomo = np.amin( np.abs( ( monoTomo, biTomo ) ) )
    vminTmp    = 20 * np.log10( max(min2D, minTomo) )
    print(f"vmin = {vminTmp:.2f}, vmax = {vmaxTmp:.2f}")

    if vmin == -1:
        vmin = vminTmp
    else:
        vmin = vmin
    if vmax == -1:
        vmax = vmaxTmp
    else:
        vmax = vmax

    # GRID 0
    toPlot = 20 * np.log10( np.abs( mono2D ) )
    im = grid0[0].imshow( np.squeeze( toPlot ), cmap='jet', extent=extentyx, vmin=vmin, vmax=vmax, origin='lower')
    grid0[0].set_ylabel("mono2D")

    toPlot = 20 * np.log10( np.abs( bi2D ) )
    im = grid0[1].imshow( np.squeeze( toPlot ), cmap='jet', extent=extentyx, vmin=vmin, vmax=vmax, origin='lower')
    grid0[1].set_ylabel("bi2D")
    grid0[1].yaxis.set_ticks_position('right')
    grid0[1].yaxis.set_label_position('right')

    # GRID 1
    toPlot = 20 * np.log10( np.abs( monoTomo ) )
    im = grid1[0].imshow( np.squeeze( toPlot ), 
        cmap='jet', extent=extentyz, vmin=vmin, vmax=vmax, origin='lower')
    grid1[0].set_ylabel("monoTomo")
    #grid[1].get_shared_x_axes().join(grid[0], grid[1])
        
    toPlot = 20 * np.log10( np.abs( biTomo ) )
    im = grid1[1].imshow( np.squeeze( toPlot ), 
        cmap='jet', extent=extentyz, vmin=vmin, vmax=vmax, origin='lower',)
    grid1[1].set_ylabel("biTomo")
    grid1[1].yaxis.set_ticks_position('right')
    grid1[1].yaxis.set_label_position('right')

    grid0.cbar_axes[0].colorbar(im)
    grid1.cbar_axes[0].colorbar(im)

    for ax in grid0:
        ax.yaxis.set_major_locator(MaxNLocator(prune='both', integer=True))

    for ax in grid1:
        ax.yaxis.set_major_locator(MaxNLocator(prune='both', integer=True))
    
    return fig, grid0, grid1
