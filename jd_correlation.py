import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import medfilt
from scipy.signal import correlate
from pylab import *
from matplotlib import rc

rc('mathtext', default='regular')


def get_xy_data(filename=None):
    """
    Utility function to retrieve x,y vectors form an ASCII data fiqle.
    It return the first column as x, the second as y.
    Usage:
        x,y = get_xy_data("mydata.txt")
    """
    data = np.loadtxt(filename, usecols=(0,1))
    return data[:,0], data[:,1]

def new_vLSR(l,v_lsr):
    Uo_IAU = 10.27            # km/s precessed to J2000
    Vo_IAU = 15.32
    Wo_IAU =  7.74
# Modern Uo, Vo, Wo values from Reid et al (2014)
    Uo = 10.50                # km/s
    Vo =  14.4
    Wo =  8.9

# Assume b=0, so
    cos_b = 1.0
    sin_b = 0.0
    v_newlsr=zeros(len(v_lsr))
    for i in range(len(v_lsr)):
        cos_l = cos(l[i]*pi/180.)  
        sin_l = sin(l[i]*pi/180.) 
        v_helio = v_lsr[i] - (Vo_IAU*sin_l + Uo_IAU*cos_l)*cos_b -  Wo_IAU*sin_b
        v_newlsr[i] = v_helio + (Vo*sin_l + Uo*cos_l)*cos_b +  Wo*sin_b
    return v_newlsr 


plt.ion()
plt.clf()
r0 = 8.5
theta0 = 220.0
l1,vel1 = get_xy_data("Q1v.dat")
l2,vel2 = get_xy_data("Q4v.dat")
v1_temp = new_vLSR(l1,vel1)    # Apply the VLSR correction
v2_temp = new_vLSR(l2,vel2)
v1 = medfilt(v1_temp,kernel_size=5)
v2 = medfilt(v2_temp,kernel_size=5)   # Do the median filter that was done in McG07
ndat1 = len(v1)
ndat2 = len(v2)
sinl1 = abs(np.sin(l1*np.pi/180))
wt1= np.array(sinl1)
sinl2 = abs(np.sin(l2*np.pi/180))
wt2= np.array(sinl2)
rot1 = np.array(abs(v1) + theta0*sinl1)
rg1 = np.array(r0*sinl1)
rot2 = np.array(abs(v2) + theta0*sinl2)
rg2 = np.array(r0*sinl2)
#
# Create arrays with both datasets
#    print "Creating array of length",ndat1+ndat2;
rg = np.append(rg1,rg2)
rot = np.append(rot1,rot2)
wt = np.append(wt1,wt2)
x = rg/r0
y = rot/theta0
weights = wt/np.mean(y)
# Fit the array (currently unweighted so dominated by dense sampling at high R)
pars = np.polyfit(x,y,1)
#    print pars
yfit = np.polyval(pars,x)
rotfit = yfit * theta0
#
# Create arrays of the differences from fit
fit1 = np.polyval(pars,rg1/r0)
fit2 = np.polyval(pars,rg2/r0)
diff1 = rot1 - fit1*theta0
diff2 = rot2 - fit2*theta0
#  the 1st quadrant data are in the opposite order now, flip both arrays
diff0 = diff1
diff1=diff0[::-1]
#
#   diff1 and diff2 are the residual velocities after subtracting the best fits
#
#  now go back to longitude for equal step in x = linear dist along circle=locus of scp
#    distance along that circle is r l (l in radians).  Note factor of 2 between
#    longitude angle and central angle on the locus of scp circle.
#
#  we'll just take the range x=3 to x=9.5 kpc (scaled by r0) defined by x01,x02
#
x0=abs(l1*r0*np.pi/180.)
#  flip this array as we flipped diff1 above
x1=x0[::-1]
x2=abs((360.-l2)*r0*np.pi/180.)
x01=np.linspace(3.,9.5,1300)
x02=np.linspace(3.,9.5,1300)
#
#  Alternative:  equal steps in radius distance, i.e. equal steps in sin(l)
#
x0=abs(r0*np.sin(l1*np.pi/180.))
#  flip this array as we flipped diff1 above
x1=x0[::-1]
x2=abs(r0*np.sin((360.-l2)*np.pi/180.))
x01=np.linspace(3.5,7.5,800)
x02=np.linspace(3.5,7.5,800)
#
#  Note the radial range is 3.5 to 7.5 kpc, with steps of 5 pc as before
#  end alternative choice
#
y01=np.interp(x01,x1,diff1)
y02=np.interp(x02,x2,diff2)
#
#  note the new step size is Delta-x = 6500/1300 pc = 5 pc
#
figure(1)
plt.plot(x1,diff1,'r+',label="QI")
plt.plot(x2,diff2,'b+',label="QIV")
plt.plot(x01,y01,'r-',label="QI interp")
plt.plot(x02,y02,'b-',label="QIV interp")
#    plt.plot(x01,y01,'ro',label="QI interp")
#    plt.plot(x02,y02,'bo',label="QIV interp")
plt.xlabel(r"$s \,(kpc)$")
plt.ylabel(r"$\Delta \Theta \, \, (km\, s^{-1})$")
#    plt.legend(loc=2)
boxx=[3.,3.,9.5,9.5,3.]
boxy=[-14.5,+14.5,+14.5,-14.5,-14.5]
plt.plot(boxx,boxy,'k-')
plt.show()
#
#-------------------------- here begins the correlation analysis
# --  z is the ccf, z1 is the acf of Q1, z2 is the acf of Q4
#
z=correlate(y01,y02)/(len(y01)*np.std(y01)*np.std(y02))
z1=correlate(y01,y01)/(len(y01)*np.var(y01))
z2=correlate(y02,y02)/(len(y01)*np.var(y02))

z1b=z1/z1[len(z1)/2]
z2b=z2/z2[len(z2)/2]
z1=z1b
z2=z2b

plt.clf()
plt.plot(x01,y01,'r.')
plt.plot(x02,y02,'b.')

#-------------------------- and now the plotting

# just for a test:
##    z=jdccf_0pad(y02,y01)/(np.std(y01)*np.std(y02))
##    z1=jdccf_0pad(y01,y01)/(np.var(y01))
##    z2=jdccf_0pad(y02,y02)/(np.var(y02))
# test done
#    print len(z)
zx1=np.linspace(1,len(z),len(z))
zx2=(zx1-(len(z)/2.))*.005
#plt.figure(2)
plt.clf()
plt.plot(zx2,z,'g-',label="ccf")
plt.plot(zx2,z1,'r-',label="acf Q1")
plt.plot(zx2,z2,'b-',label="acf Q4")
plt.legend(loc=1)
##    plt.xlabel(r"$R \,(kpc)$")
##    plt.ylabel(r"$\Delta\Theta (km\, s^{-1})$")
plt.xlabel(r"$\Delta x \, \, \, (kpc)$")
#  use the above for distance along lscp circle computation
#  use the below for radial distance computation
plt.xlabel(r"$\Delta r_G \, \, \, \mathrm{(kpc)}$")
plt.ylabel("Normalised Correlation Coefficient")
plt.minorticks_on()
plt.vlines(0.,-1.5,1.5,linestyles='dotted',linewidth=1.2)
plt.ylim(-0.5,1.2)

#    plt.xlim(-3.5,3.5)
#
#  the y axis is normalized to 1 for 100% correlation
#  the x axis is in lag steps of 7.2 pc
#
#-------------------------- all done


#     keep this for later, in case plots of the raw residuals vs. R are needed
#    plt.plot(rg1,diff1,'r+',label="QI")
#    plt.plot(rg2,diff2,'b+',label="QIV")
#    plt.xlabel(r"$R \,(kpc)$")
#    plt.ylabel(r"$\Delta\Theta (km\, s^{-1})$")
#    plt.xlim(3.0,8.0)
#    plt.ylim(-15.,15.)
#    print np.mean(diff1),np.mean(diff2)
#    print np.median(diff1),np.median(diff2)
#    print np.std(diff1),np.std(diff2)
#    plt.legend(loc=2)
#    plt.show()

