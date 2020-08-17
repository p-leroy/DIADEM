import configparser, os

import numpy as np

# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 09:34:38 2020

@author: pleroy
"""

class Conf(object):
    def __init__(self):
        self.nb_elev = None
        self.nb_freq = None
        self.nb_ssb = None
        self.el0 = None
        self.az0 = None
        self.half_range = None
        self.elevation = None

class ConfigSer2(object):
    def __init__(self, filename):
        self.filename = filename
        self.ref = Conf()
        self.sup = Conf()
        self.load()
        
    def load(self):
        config = configparser.ConfigParser()
        config.read(self.filename)
        
        # ['reference']
        self.ref.nb_elev = int(config['reference']['nb_elev'])
        self.ref.nb_freq = int(config['reference']['nb_freq'])
        self.ref.nb_ssb = int(config['reference']['nb_ssb'])
        self.ref.el0 = float(config['reference']['el0'])
        self.ref.az0 = float(config['reference']['az0'])
        self.ref.half_range = float(config['reference']['half_range'])
        self.ref.elevation = np.linspace(
            self.ref.el0 - self.ref.half_range, 
            self.ref.el0 + self.ref.half_range, 
            self.ref.nb_elev)
        
         # ['support']
        self.sup.nb_elev = int(config['support']['nb_elev'])
        self.sup.nb_freq = int(config['support']['nb_freq'])
        self.sup.nb_ssb = int(config['support']['nb_ssb'])
        self.sup.el0 = float(config['support']['el0'])
        self.sup.az0 = float(config['support']['az0'])
        self.sup.half_range = float(config['support']['half_range'])
        self.sup.elevation = np.linspace(
            self.sup.el0 - self.sup.half_range, 
            self.sup.el0 + self.sup.half_range, 
            self.sup.nb_elev)
        
        # ['data']
        self.dic_dir = {key: eval(config['data'][key]) for key in config['data']}
        
        # ['path']
        self.base_path = os.path.dirname(self.filename)
        self.out_path = os.path.join(self.base_path, config['path']['out_path'])
        
        # ['processing']
        self.gateWidth = float(config['processing']['gateWidth'])
        self.gateCenter = float(config['processing']['gateCenter'])
        self.centerFreq = float(config['processing']['centerFreq'])
        self.bandWidth = float(config['processing']['bandWidth'])
        self.idxElevation = float(config['processing']['idxElevation'])
        self.vmin = float(config['processing']['vmin'])
        self.vmax = float(config['processing']['vmax'])

    def printConfig(self):
        print(self.filename)
        print(f"base_path {self.base_path}")
        print(f"out_path {self.out_path}")
        print(f"gateWidth {self.gateWidth}")
        print(f"gateCenter {self.gateCenter}")
        print(f"centerFreq {self.centerFreq}")
        print(f"bandWidth {self.bandWidth}")
        print(f"idxElevation {self.idxElevation}")
        print(f"vmin {self.vmin}")
        print(f"vmax {self.vmax}")
        print(f"\nSUPPORT")
        print(f"nb_elev {self.sup.nb_elev}")
        print(f"nb_freq {self.sup.nb_freq}")
        print(f"nb_ssb {self.sup.nb_ssb}")
        print(f"el0 {self.sup.el0}")
        print(f"az0 {self.sup.az0}")
        print(f"\nREFERENCE")
        print(f"nb_elev {self.ref.nb_elev}")
        print(f"nb_freq {self.ref.nb_freq}")
        print(f"nb_ssb {self.ref.nb_ssb}")
        print(f"el0 {self.ref.el0}")
        print(f"az0 {self.ref.az0}")

