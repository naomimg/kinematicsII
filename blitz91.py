import numpy as np
from pylab import *
from scipy.signal import medfilt
from matplotlib import rc
from scipy.signal import resample
from scipy.interpolate import interp1d

rc('mathtext', default='regular')

# Plot the Blitz & Spergel 1991 function

def Theta(l,theta_b,switch=True):
	epsilon = 0.02
	q = 3.0
	R0 = 8000.
	vc = 220.0
# switch let's us choose the positive or negative variation
	if switch:
		little_theta = math.pi/2. + np.radians(l) + np.radians(theta_b)
	else:
		little_theta = math.pi/2. - np.radians(l) + np.radians(theta_b)
	sinl  = np.sin(np.radians(l))

	Thetasq = R0 * sinl * ((1./(R0*sinl)) + epsilon/2. * ((3* np.cos(little_theta)**2 - 1) * (2-q) * sinl**(1-q) / R0))
	Theta = sqrt(Thetasq)* vc
	return Theta

# Overplot Blitz \& Spergel (1991)
Pi0=14.0
l_mod = arange(0.,90.)
coslterm = 2.* Pi0 * np.cos(np.radians(l_mod))
sinl = np.sin(np.radians(l_mod))
for theta_b in (30.,40.,50.,60.):
	a = Theta(l_mod,theta_b,switch=True)
	b = Theta(l_mod,theta_b,switch=False)
	y = a-b
	deltav = coslterm - y
	plot(sinl,y)
#plot(l_mod,coslterm,'k')
#plot(l_mod,Theta(l_mod,theta_b,vc,switch=True),'g')
#plot(l_mod,Theta(l_mod,theta_b,vc,switch=False),'r')
plot(l_mod,y,'purple')


xlim(20,65.)
ylim(-20.,25.)
minorticks_on()
pylab.text(40,21.4,'$\Pi_0=14$')
pylab.text(25,18.3,'$\Pi_0=10$')
pylab.text(23,13,'$\Pi_0=7$')
savefig("vdiff_vs_l.pdf")
