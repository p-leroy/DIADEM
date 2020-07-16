import configparser

import numpy as np

# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 09:34:38 2020

@author: pleroy
"""

class ConfigSer(object):
    def __init__(self, filename):
        self.filename = filename
        self.load()
        
    def load(self):
        config = configparser.ConfigParser()
        config.read(self.filename)
        
        # ['parameters']
        self.nb_elev = int(config['parameters']['nb_elev'])
        self.nb_freq = int(config['parameters']['nb_freq'])
        self.nb_ssb = int(config['parameters']['nb_ssb'])
        self.el0 = float(config['parameters']['el0'])
        self.az0 = float(config['parameters']['az0'])
        self.half_range = float(config['parameters']['half_range'])
        self.elevation = np.linspace(
            self.el0 - self.half_range, 
            self.el0 + self.half_range, 
            self.nb_elev)
        self.base_path = config['parameters']['base_path']
        self.out_path = config['parameters']['out_path']
        
        # ['data']
        self.dic_dir = {key: eval(config['data'][key]) for key in config['data']}
        
        # ['gating']
        self.gateWidth = float(config['gating']['gateWidth'])
        self.gateCenter = float(config['gating']['gateCenter'])

    def printConfig(self):
        print(self.filename)
        print(f"nb_elev {self.nb_elev}")
        print(f"nb_freq {self.nb_freq}")
        print(f"nb_ssb {self.nb_ssb}")
        print(f"el0 {self.el0}")
        print(f"az0 {self.az0}")
        print(f"base_path {self.base_path}")
        print(f"out_path {self.out_path}")
        print(f"gateWidth {self.gateWidth}")
        print(f"gateCenter {self.gateCenter}")
