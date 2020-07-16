import numpy as np

# -*- coding: utf-8 -*-
"""
Éditeur de Spyder

Ceci est un script temporaire.
"""

class Probe(object):
    def __init__(self, filename):
        self.filename = filename
        self.load()
        
    def load(self):
        tmp = np.loadtxt(self.filename, skiprows=5)
        self.S11 = (10 ** (0.05 * tmp[:,1])) * np.exp(1j * np.pi / 180 * tmp[:,2])
        self.freq = tmp[:,0] / 1e9
        self.Nfreq = self.freq.size
        self.d = np.linspace( 0, 300 / (self.freq[1] - self.freq[0]), self.freq.size) * 0.5
        self.data = tmp
