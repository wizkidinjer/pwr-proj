# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from scipy.special import jv
import matplotlib.pyplot as plt
import numpy as np


''''REACTOR CONFIGURATION'''

EperF=200 * 1.6022e-13 #joules per fission





class pwr:
    fuelrho=8e20 #at/cm3
    sigf= 585e-24 #cm2
    P= 500e6
    f0=5e13
    

    def __init__(self, R, Z,):
        """
        R       - extrapolated radius (cm)
        Z       - extrapolated height (cm)
        Power   - reactor power (W)
        fuelrho - fuel number density (at/cm3)
        sigf    - fission cross section (cm2)
        """
        self.R = R
        self.Z = Z

    def get_dims(self, R, Z,):
        


    def flux(self, Z,R):
        h=np.linspace(-L/2,L/2, int(1e4))
        r=np.linspace(0,R, int(1e3))
        
        #constant properties
        
        f= f0 * jv(0,2.405*r/R) * np.cos(np.pi * h/L)
        powerdens = sigf * EperF f
        
        Power=