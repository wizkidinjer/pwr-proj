# -*- coding: utf-8 -*-


from scipy.special import jv
import matplotlib.pyplot as plt
import numpy as np
from numpy import pi
import scipy as si


''''REACTOR CONFIGURATION'''

EperF=200 * 1.6022e-13 #joules per fission





class pwr:
    fuelrho = 8e20        # at/cm3
    sigf    = 585e-24     # cm2
    f0      = 5e13        # peak flux (n/cm2/s)
    EperF   = 3.2e-11     # ~200 MeV in joules

    def __init__(self, R, H):
        self.R = R   # extrapolated radius (cm)
        self.H = H   # extrapolated height (cm)

    def flux(self, r, h):
        """neutron flux at a single (r, h) point"""
        return self.f0 * jv(0, 2.405 * r / self.R) * np.cos(pi * h / self.H)

    def power(self):

        integrand = lambda h, r: self.fuelrho * self.sigf * EperF * self.f0 * jv(0, 2.405 * r / self.R) * np.cos(pi * h / self.H)

        P, err = si.integrate.dblquad(
            integrand,
            0, self.R,-self.H / 2,self.H / 2  )
        return P, err
    

core=pwr(R=150,H=400)
print(core.power())


r=np.linspace(0,core.R,1000)
h=np.linspace(-core.H/2,core.H/2,1000)
rgrid,hgrid=np.meshgrid(r,h)
f=core.flux(rgrid,hgrid)

from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(rgrid, hgrid, f, cmap='inferno',)
ax.set_xlabel('r (cm)')
ax.set_ylabel('h (cm)')
ax.set_zlabel('Flux')
ax.set_title('Flux Distribution')
plt.show()


