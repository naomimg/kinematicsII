import galpy.potential
import numpy as np
import pylab
import math

from galpy.potential import MWPotential2014
from galpy.potential import evaluatePotentials as evalPot
from galpy.orbit import Orbit

from optparse import OptionParser


from galpy.potential import SteadyLogSpiralPotential, \
    LogarithmicHaloPotential, TransientLogSpiralPotential, \
    DehnenBarPotential, PowerSphericalPotential, \
    MovingObjectPotential, lindbladR, EllipticalDiskPotential
from galpy.df_src.evolveddiskdf import evolveddiskdf
from galpy.util import save_pickles, bovy_plot, bovy_coords, multi
from galpy.orbit import Orbit
from galpy.df import dehnendf, shudf

from galpy.util import config

baralpha=0.01  #Bar strength alpha parameter
bar_tform=-4.  # Time for Bar to form in bar periouds
bar_tsteady=2.0  # Bar steady time in bar periods
beta=0.0 # Rotation curve power law index
bar_olr=.9  # Radius of OLR
bar_angle=25.0 # Bar angle
config_setro=8.
axip= LogarithmicHaloPotential(normalize=1.,q=0.95)  # Set up first potential
pot= [axip]
barp= DehnenBarPotential(alpha=baralpha,
                                 tform=bar_tform,
                                 tsteady=bar_tsteady,
                                 beta=beta,
                                 rolr=bar_olr,
                                 barphi=bar_angle/180.*np.pi,
                                 ro=8.0,
                                 vo=220.)
# Add Bar

pot.append(barp)
to= barp.tform()


# Evolve the disk
rd=1./3.
rs=2.0
so = 31.4/220.
dfc= dehnendf(beta=beta,correct=True,niter=20,profileParams=(rd,rs,so),savedir='./')

edf= evolveddiskdf(dfc,pot,to=to)

#vt_n = galpy.potential.vterm(EllipticalDiskPotential,l_deg_n) * 220.0

# Plot the surfacemass for the distribution function
Rs= np.linspace(0.01,5.,151)
out= [dfc.surfacemass(r) for r in Rs]
pylab.plot(Rs, out)

# Setup a grid
xgrid= np.linspace(0.,2.*math.pi*(1.-1./26./2.),
                                   2.*26.)
ygrid= np.linspace(0.1,2.0,26.)

nx= len(xgrid)
ny= len(ygrid)
print nx,ny
ii, jj= 0, 0
nt=1
surfmass= np.zeros((nx,ny,1))
meanvr= np.zeros((nx,ny,1))
meanvt= np.zeros((nx,ny,1))
grids= []
evalts= np.linspace(0.,to,nt)
nsigma=1
R, phi= ygrid[jj], xgrid[ii]
while  ii < nx:
    #Calculate surfmass etc.
    while jj < ny:
        smass, grid= edf.vmomentsurfacemass(R,0,0,grid=True,phi=phi,
                                        returnGrid=True,t=evalts,
                                        gridpoints=26,
                                        nsigma=nsigma,
                                        print_progress=True)
        surfmass= smass
        meanvr[ii,jj,:]= edf.meanvR(R,phi=phi,grid=grid,t=evalts,
                       surfacemass=surfmass,
                       nsigma=nsigma)
        meanvt[ii,jj,:]= edf.meanvT(R,phi=phi,grid=grid,t=evalts,
                       surfacemass=surfmass,
                       nsigma=nsigma)


# Pickle the output
savefilename = 'output.pickle'
savefile= open(savefilename,'wb')
pickle.dump(surfmass,savefile)
pickle.dump(meanvr,savefile)
pickle.dump(meanvt,savefile)

