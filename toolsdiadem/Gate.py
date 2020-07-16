import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import kaiser

# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 09:12:35 2020

@author: pleroy
"""


class Gate(object):
    def __init__(self, d, width=500, center=6650):
        self.d = d
        self.N = d.size
        self.width = width
        self.center = center
        
        self.buildGate()
        self.buildGateKB()
        self.buildGateKB2()
        self.buildGateLLC()

    def buildGate(self):
        self.gate = np.zeros(self.N)
        self.width = 500
        self.size = np.where(np.abs(self.d - self.center) < self.width)[0].shape[0]
        self.gate[np.where(np.abs(self.d - self.center) < self.width)] = 1.0
        self.fft = np.fft.fft(self.gate)
        print(f"gateSize {self.size}")

    def buildGateKB(self):
        # scipy.signal.kaiser(M, beta, sym=True)
        # **SYM** When True (default), generates a symmetric window, for use in filter design.
        # When False, generates a periodic window, for use in spectral analysis
        beta = 6
        myWindow = np.zeros(self.N)
        myWindow[0] = 1
        myWindow[1:] = np.fft.fftshift(kaiser(self.N-1, beta, sym=True))
        modMyWindow = 20 * np.log10(np.abs(myWindow))
        self.fftKB = np.fft.fft(self.gate) * myWindow
        self.gateKB = np.fft.ifft(self.fftKB)

    def buildGateKB2(self):
        NKaiser = 200
        littleKaiser = kaiser(2 * NKaiser, 6, sym=True)
        myWindowLittleKaiser = np.zeros(self.N)
        myWindowLittleKaiser[0] = 1
        myWindowLittleKaiser[1:NKaiser + 1] = littleKaiser[NKaiser:]
        myWindowLittleKaiser[self.N - NKaiser:] = littleKaiser[0:NKaiser]
        fftMyGateLittleKaiser = self.fft * myWindowLittleKaiser
        self.gateKB2 = np.fft.ifft(fftMyGateLittleKaiser)

    def buildGateLLC(self, id="llc"):
        gate = 0 * self.d + 1e-6
        gate[np.where(np.abs(self.d - self.center) < self.width)] = 1.0
        
        d_centre = np.sum(gate * self.d) / np.sum(gate)
        
        gate_a0 = (1e-6 + ((self.d - d_centre) / d_centre) ** 4)
        gate_a = 1 / gate_a0
        
        offset = np.amax(gate_a * (1 - gate))
        gate_b0 = gate_a / offset - 1.0
        gate_b  = gate_b0 * (1 - gate) + 1.0
        
        tmp = np.where(gate == 1)[0]
        y1 = self.d * 0.0
        y1[0 : tmp[0]] = self.d[0:tmp[0]]
        y2 = self.d * 0.0
        y2[tmp[-1]:] = self.d[tmp[-1]:]
        
        y1 = y1 / self.d[tmp[0]]
        y2[tmp[-1]:] = (y2[tmp[-1]:] - np.amax(self.d)) / (self.d[tmp[-1]] - np.amax(self.d))
        
        y1[tmp[0]] = 0.0
        y2[tmp[-1]] = 0.0
        
        gate_c = gate_b * (y1 + y2 + gate)
        self.gateLLC = gate_c
        
        return gate
    
    def plotGateLLC(self):
        fig, ax = plt.subplots(1,1)
        ax.plot(d, gate, label='gate')
        ax.legend()
        ax.grid()
        fig.savefig("gate.png")
