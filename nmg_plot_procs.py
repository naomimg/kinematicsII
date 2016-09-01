#!/usr/bin/env python
import numpy as np
from pylab import *
from scipy.signal import medfilt
from matplotlib import rc
from scipy.signal import resample
from scipy.interpolate import interp1d

rc('mathtext', default='regular')

def get_xy_data(filename=None):
    """
    Utility function to retrieve x,y vectors from an ASCII data file.
    It returns the first column as x, the second as y.
    Usage:
        x,y = get_xy_data("mydata.txt")
    """
    import numpy as np
    data = np.loadtxt(filename, usecols=(0,1))
    return data[:,0], data[:,1]

def plot_dens(file=None):
    import matplotlib.pyplot as plt
    plt.ion()
    plt.clf()
    if file is None:
        raise ArgumentError("please give file name to plot")
    x, y = get_xy_data(file)
    plt.plot(x,y)
    xlbl = r"$R_g (\rm{kpc})$"
    plt.xlabel(xlbl)
    ylbl = r"$n_{HI} \ (\rm{cm^{-3}})$"
    plt.ylabel(ylbl)
    plt.ylim(-0.05,4)
    plt.xlim(0,25)
    plt.show()


def plot_vels(velsfile=None,plot_label=None):
    import matplotlib.pyplot as plt
    plt.ion()
    plt.clf()
    if velsfile is None:
        raise ArgumentError("Please give a file name to plot")
    l, v = get_xy_data(velsfile)
    plt.plot(l,v,'k+',label=plot_label,alpha=1)
#    xlim(340.,288.)
    plt.legend()
    plt.xlabel(r'$Galactic\, Longitude\, (deg)$')
    plt.ylabel(r'$V_t\, (km\, s^{-1})$')
    plt.show()

def vsinl_plots(velsfile1=None,velsfile2=None):
    import matplotlib.pyplot as plt
    plt.ion()
    plt.clf()
    import numpy as np
    from scipy.signal import medfilt
    if velsfile1 is None:
        raise ArgumentError("Please give input file")
    l1,v1 = get_xy_data(velsfile1)
    l2,v2 = get_xy_data(velsfile2)
    sinl1 = abs(np.sin(l1*np.pi/180))
    sinl2 = abs(np.sin(l2*np.pi/180))
    v2_med = medfilt(v2,kernel_size=5)
    v2 = abs(v2_med)
    plt.plot(sinl1,v1,'b+',label="QI")
    plt.plot(sinl2,v2,'r+',label="QIV")
    plt.legend()
    plt.xlabel(r'$sin(l)$')
    plt.ylabel(r'$|V_t|\, (km\, s^{-1})$')
    # Fit the arrays
    x = sinl1
    y = v1
    pars = np.polyfit(x,y,1)
    print pars
    yfit = np.polyval(pars,x)
    plt.plot(sinl1,yfit,'b--')
    x = sinl2
    y = v2
    pars = np.polyfit(x,y,1)
    print pars
    yfit = np.polyval(pars,x)
    plt.plot(sinl2,yfit,'r--')
    plt.ylim(15.,140.)
    plt.show()


def plot_vels2(velsfile1=None,velsfile2=None):
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator
    plt.ion()
    plt.clf()
    import numpy as np
    from scipy.signal import medfilt
    if velsfile1 is None:
        raise ArgumentError("Please give a file name to plot")
# Put the minor tick marks on the plot
    xminorLocator=MultipleLocator(2)
    yminorLocator=MultipleLocator(5)
    ax1=plt.subplot(111)
    l1,v1 = get_xy_data(velsfile1)
    l2,v2 = get_xy_data(velsfile2)
    v2 = medfilt(v2,kernel_size=3)
    v1 =  medfilt(v1,kernel_size=3)
#    l2=abs(360.0-l2)
    v2 = -v2
    plt.plot(l1,v1,'r+',label="QI",alpha=1)
    ax1.xaxis.tick_bottom()   
    plt.xlabel(r'$Galactic\, Longitude\, (deg)$')
    plt.ylabel(r'$|V_t|\, (km\, s^{-1})$')
    plt.xlim(15.,70.)
    from matplotlib import rc
    rc('mathtext', default='regular')
    ax1.xaxis.set_minor_locator(xminorLocator) 
    ax1.yaxis.set_minor_locator(yminorLocator) 
    plt.minorticks_on()
    # Fit the array (currently unweighted so dominated by dense sampling at high R)
    x = abs(np.sin(l1*np.pi/180))
    y = v1
    pars = np.polyfit(x,y,1)
    print pars
    yfit = np.polyval(pars,x)
    plt.plot(l1,yfit,'r--')
    ax2= plt.twiny()
    ax2.xaxis.tick_top()   
    plt.plot(l2,v2,'b+',label="QIV",alpha=1)
    plt.xlim(345.,290.)
    ax2.xaxis.set_minor_locator(xminorLocator) 
    ax2.yaxis.set_minor_locator(yminorLocator) 
        # Fit the array (currently unweighted so dominated by dense sampling at high R)
    x = abs(np.sin(l2*np.pi/180))
    y = v2
    pars = np.polyfit(x,y,1)
    print pars
    yfit = np.polyval(pars,x)
    plt.plot(l2,yfit,'b--')
    plt.show()
    plt.figure(2)
    x1 = abs(np.sin(l1*np.pi/180))
    resid1 = np.polyval(pars,x1) - v1
    x2 = abs(np.sin(l2*np.pi/180))
    l2a = abs(l2-360.)
    resid2 = np.polyval(pars,x2) - v2
    plt.plot(l1,resid1,'r+')
    plt.plot(l2a,resid2,'b+')
    plt.xlabel(r"$|l|\, (deg)$")
    plt.ylabel(r"$\Delta V_t\, (km\,s^{-1})$")
    plt.show()

def rot_plots(velsfile1=None,velsfile2=None):
    import numpy as np
    import matplotlib.pyplot as plt
    plt.ion()
    plt.clf()
    if velsfile1 is None:
        raise ArgumentError("Please give input file")
    from scipy.signal import medfilt
    r0 = 8.5
    theta0 = 220.0
#    r0=8.34    # Reid et al (2014) values
#    theta0=240.0
    l1,vel1 = get_xy_data(velsfile1)
    l2,vel2 = get_xy_data(velsfile2)
    v1 = medfilt(vel1,kernel_size=5)
    v2 = medfilt(vel2,kernel_size=5)   # Do the median filter that was done in McG07
    sinl1 = abs(np.sin(l1*np.pi/180))
    sinl2 = abs(np.sin(l2*np.pi/180))
    rot1 = abs(vel1) + theta0*sinl1
    rg1 = r0*sinl1
    rot2 = abs(v2) + theta0*sinl2
    rg2 = r0*sinl2
    plt.plot(rg1,rot1,'r+',label="QI")
    plt.plot(rg2,rot2,'b+',label="QIV")
    plt.minorticks_on()
    plt.legend()
    plt.xlabel(r"$R\, (kpc)$")
    plt.ylabel(r"$\Theta\, (km\,s^{-1})$")
    print 'Number of 1st quad points: ',len(rot1) 
    print 'Number of 1st quad points: ',len(rot2)
    plt.xlim(3.0,8.0)
    plt.ylim(190.,260.)
    plt.show()

def rot_fit(velsfile1=None,velsfile2=None):
    import numpy as np
    import matplotlib.pyplot as plt
    plt.ion()
    plt.clf()
    import numpy as np
    if velsfile1 is None:
        raise ArgumentError("Please give input file")
    from scipy.signal import medfilt
    r0 = 8.5
    theta0 = 220.0
    l1,vel1 = get_xy_data(velsfile1)
    l2,vel2 = get_xy_data(velsfile2)
    v1 = new_vLSR(l1,vel1)    # Apply the VLSR correction
    v2_temp = new_vLSR(l2,vel2)
    #v1 = medfilt(vel1,kernel_size=5)
    v2 = medfilt(v2_temp,kernel_size=5)   # Do the median filter that was done in McG07
    #v1 = medfilt(vel1,kernel_size=5)
    #v1 = vel1
    #v2 = medfilt(vel2,kernel_size=5)# Do the median filter that was done in McG07
    ndat1=len(v1)
    ndat2 = len(v2)
    sinl1 = abs(np.sin(l1*np.pi/180))
    wt1= array(sinl1)
    sinl2 = abs(np.sin(l2*np.pi/180))
    wt2=array(sinl2)
    rot1 = array(abs(v1) + theta0*sinl1)
    rg1 = array(r0*sinl1)
    rot2 = array(abs(v2) + theta0*sinl2)
    rg2 = array(r0*sinl2)
#
# Create arrays with both datasets
    print "Creating array of length",ndat1+ndat2;
    rg = append(rg1,rg2)
    rot = append(rot1,rot2)
    wt = append(wt1,wt2)
    x = rg/r0
    y = rot/theta0
    weights = wt/mean(y)
# Fit the array (currently unweighted so dominated by dense sampling at high R)
    pars = polyfit(x,y,1)
    print pars
    yfit = polyval(pars,x)
    rotfit = yfit * theta0
#
# Create arrays of the differences from fit
    fit1 = polyval(pars,rg1/r0)
    fit2 = polyval(pars,rg2/r0)
    diff1 = rot1 - fit1*theta0
    diff2 = rot2 - fit2*theta0
#
# Show the B&B93 fit 
    p=[1.00767, 0.0394, 0.000712]
    bbr = arange(3.0,8.0,0.01)
    bbx = bbr/r0
    bbfit = (p[0]*bbx**p[1] + p[2])*theta0
#
# Clemens 1985 curve
    cra=arange(3.825,13.6,0.01)
    crb=arange(0.765,3.825,0.01)
    crota=np.zeros(len(cra),float)
    pa=[-2342.65,2507.60,-1024.06,224.5627,-28.40800,2.0697,-0.080508,0.00129]
    pa.reverse()
    crota=polyval(pa,cra)
#    polya=poly1d(pa)
#    crota=polya(cra)
    print max(crota),min(crota)     
    pb=[325.09,-248.14,231.87, -110.73,25.07,-2.11]
    pb.reverse()
    crotb=polyval(pb,crb)
    crot = append(crotb,crota)
    cr = append(crb,cra)
    crot = crot *226.0/220.0
#
# Show the Reid et al (2014) curve
    rrg = arange(3.,8.0,0.1)
    rr = theta0 + 0.2 * (rrg - r0)
# Plot the points and fits
    plt.plot(rg1,rot1,"r+",label="QI")
    plt.plot(rg2,rot2,"b+",label="QIV")
    plt.plot(rg,rotfit,label="Joint fit",marker="None",color="black",linestyle="-")
    plt.plot(bbr,bbfit,linestyle=":",color="black",label="BB93 fit")
    #plot(cr,crot,linestyle="-.",color="black",label="Clemens (1985)")
    plt.plot(rrg,rr,linestyle="-.",color="black",label="Reid14")
    plt.legend(loc=2)
    plt.minorticks_on()
    plt.xlabel(r"$R\, (kpc)}$")
    plt.ylabel(r"$\Theta\, (km\, s^{-1})$")
    print 'Number of 1st quad points: ',len(rot1) 
    print 'Number of 4th quad points: ',len(rot2)
    plt.xlim(3.0,8.0)
    plt.ylim(180.0,267.0)
    figure(2)
    plt.plot(rg1,diff1,'r+',label="QI")
    plt.plot(rg2,diff2,'b+',label="QIV")
    plt.xlabel(r"$R \,(kpc)$")
    plt.ylabel(r"$\Delta\Theta (km\, s^{-1})$")
    plt.xlim(3.0,8.0)
    plt.ylim(-15.,15.)
    print mean(diff1),mean(diff2)
    print median(diff1),median(diff2)
    print std(diff1),std(diff2)
    legend(loc=2)
    plt.show()
    return

def new_vLSR(l,v_lsr):
    import numpy as np
    Uo_IAU = 10.27            # km/s precessed to J2000
    Vo_IAU = 15.32
    Wo_IAU =  7.74
# Modern Uo, Vo, Wo values from Reid et al (2014)
    Uo = 10.00                # km/s
    Vo =  5.25
    Wo =  7.17
#
# Assume b=0, so
    cos_b = 1.0
    sin_b = 0.0
    v_newlsr=np.zeros(len(v_lsr))
    for i in range(len(v_lsr)):
        cos_l = np.cos(l[i]*np.pi/180.)  
        sin_l = np.sin(l[i]*np.pi/180.) 
        v_helio = v_lsr[i] - (Vo_IAU*sin_l + Uo_IAU*cos_l)*cos_b -  Wo_IAU*sin_b
        v_newlsr[i] = v_helio + (Vo*sin_l + Uo*cos_l)*cos_b +  Wo*sin_b
    return v_newlsr 

def rot_corr(velsfile1=None,velsfile2=None):
    import matplotlib.pyplot as plt
    plt.ion()
    plt.clf()
    import numpy as np
    if velsfile1 is None:
        raise ArgumentError("Please give input file")
    from scipy.signal import medfilt
    r0 = 8.5
    theta0 = 220.0
#    r0=8.34    # Reid et al (2014) values
#    theta0=240.0
    l1,vel1 = get_xy_data(velsfile1)
    l2,vel2 = get_xy_data(velsfile2)
    v1_temp = new_vLSR(l1,vel1)    # Apply the VLSR correction
    v2_temp = new_vLSR(l2,vel2)
    v1 = medfilt(v1_temp,kernel_size=5)
    v2 = medfilt(v2_temp,kernel_size=5)   # Do the median filter that was done in McG07
    sinl1 = abs(np.sin(l1*np.pi/180))
    sinl2 = abs(np.sin(l2*np.pi/180))
    rot1 = np.abs(v1) + theta0*sinl1
    rg1 = r0*sinl1
    rot2 = np.abs(v2) + theta0*sinl2
    rg2 = r0*sinl2
    plt.plot(rg1,rot1,'r+',label="QI")
    plt.plot(rg2,rot2,'b+',label="QIV")
    plt.minorticks_on()
    plt.legend()
    plt.xlabel(r"$R\, (kpc)$")
    plt.ylabel(r"$\Theta\, (km\,s^{-1})$")
    print 'Number of 1st quad points: ',len(rot1) 
    print 'Number of 1st quad points: ',len(rot2)
    plt.xlim(3.0,8.0)
    plt.ylim(190.,260.)
    plt.show()


def vels2_corr(velsfile1=None,velsfile2=None):
    import matplotlib.pyplot as plt
    plt.ion()
    plt.clf()
    from scipy.signal import medfilt
    if velsfile1 is None:
        raise ArgumentError("Please give a file name to plot")
# Put the minor tick marks on the plot
    xminorLocator=plt.MultipleLocator(2)
    yminorLocator=plt.MultipleLocator(5)
    ax1=plt.subplot(111)
    l1,vel1 = get_xy_data(velsfile1)
    l2,vel2 = get_xy_data(velsfile2)
#
    v1 = new_vLSR(l1,vel1)    # Apply the VLSR correction
    v2_temp = new_vLSR(l2,vel2)
    v2 = medfilt(v2_temp,kernel_size=3)
#    l2=abs(360.0-l2)
    v2 = -v2
    vel2 = -vel2
    l2 = abs(l2-360.)
    plt.plot(l1,v1,'r+',label="QI",alpha=1)
    plt.plot(l2,v2,'b+',label="QIV",alpha=1)
    plt.plot(l1,vel1,'ro',label="QI",alpha=1)
    plt.plot(l2,vel2,'bo',label="QIV",alpha=1)


def jd_rplots(velsfile1="Q1v.dat",velsfile2="Q4v.dat"):
    import numpy as np
    import matplotlib.pyplot as plt
    plt.ion()
    plt.clf()
    if velsfile1 is None:
        raise ArgumentError("Please give input file")
    from scipy.signal import medfilt
    r0 = 8.5
    theta0 = 220.0
#    r0=8.34    # Reid et al (2014) values
#    theta0=240.0
    l1,vel1 = get_xy_data(velsfile1)
    l2,vel2 = get_xy_data(velsfile2)
    v1 = medfilt(vel1,kernel_size=5)
    v2 = medfilt(vel2,kernel_size=5)   # Do the median filter that was done in McG07
    sinl1 = abs(np.sin(l1*np.pi/180))
    sinl2 = abs(np.sin(l2*np.pi/180))
    xl1 = abs(l1*np.pi/180)
    xl2 = abs((l2-360.)*np.pi/180)
    rot1 = abs(vel1) + theta0*sinl1
    rg1 = r0*sinl1
    xg1 = r0*xl1
    rot2 = abs(v2) + theta0*sinl2
    rg2 = r0*sinl2
    xg2 = r0*xl2
#  old version    plt.plot(rg1,rot1,'r+',label="QI")
#  old version    plt.plot(rg2,rot2,'b+',label="QIV")
    plt.plot(xg1,rot1,'r+',label="QI")
    plt.plot(xg2,rot2,'b+',label="QIV")
    plt.minorticks_on()
    plt.legend(loc=4)
#  old version    plt.xlabel(r"$R\, (kpc)$")
    plt.xlabel(r"$x\, (kpc)$")
    plt.ylabel(r"$\Theta\, (km\,s^{-1})$")
    print 'Number of 1st quad points: ',len(rot1) 
    print 'Number of 1st quad points: ',len(rot2)
#    plt.xlim(3.0,8.0)
#    plt.ylim(190.,260.)
    plt.show()

# ------------------------------------------
#  this is just for testing to confirm scipy.signal.correlate is working
#  it is
#
def jdccf_0pad(x,y):
  import numpy as np
  if (len(x) != len(y)):
    print "error unequal input lengths"
  if (len(x) < 2 ):
    print "error too short input lengths"
  n=len(x)-1
  xret=[]
  xnorm=[]
  x1=x
  y1=y
  for i in range(len(x)):
     z=x1*y1
     xret.insert(i,sum(z))
     xnorm.insert(i,len(x)-i)
     x2=np.roll(x1,1)
     x2[0]=0
     x1=x2
  xr1=np.divide(xret,xnorm)
  print xnorm[0:10]
  print xret[0:10]
  x1=x
  xret=[]
  xnorm=[]
  for i in range(1,len(x)):
     x2=np.roll(x1,-1)
     x2[n]=0
     x1=x2
     z=x1*y1
     xret.insert(i,sum(z))
     xnorm.insert(i,len(x)-i)
  xr2=np.divide(xret,xnorm)
  print xnorm[0:10]
  print xret[0:10]
  xr0=xr2[::-1]
  xret=np.concatenate((xr0,xr1))
  return xret

# ------------------------------------------

def jd_rot_fit(velsfile1="Q1v.dat",velsfile2="Q4v.dat"):
#
#  this version uses x as the spatial variable (LSCP length)
#
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.signal import medfilt,correlate
    plt.ion()
    plt.clf()
    if velsfile1 is None:
        raise ArgumentError("Please give input file")
    r0 = 8.5
    theta0 = 220.0
    l1,vel1 = get_xy_data(velsfile1)
    l2,vel2 = get_xy_data(velsfile2)
    v1_temp = new_vLSR(l1,vel1)    # Apply the VLSR correction
    v2_temp = new_vLSR(l2,vel2)
    v1 = medfilt(v1_temp,kernel_size=5)
    v2 = medfilt(v2_temp,kernel_size=5)   # Do the median filter that was done in McG07
    ndat1=len(v1)
    ndat2 = len(v2)
    sinl1 = abs(np.sin(l1*np.pi/180))
    wt1= np.array(sinl1)
    sinl2 = abs(np.sin(l2*np.pi/180))
    wt2=np.array(sinl2)
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
#
#   diff1 and diff2 are the residual velocities after subtracting the best fits
#
#  now go back to longitude for equal step in x = linear dist along magic circle
#  we'll just take the range x=3 to x=9.5 kpc (scaled by r0) defined by x01,x02
#
    x0=abs(l1*r0*np.pi/180.)
#  the 1st quadrant data are in the opposite order now, flip both arrays
    diff0 = diff1
    x1=x0[::-1]
    diff1=diff0[::-1]
    x2=abs((360.-l2)*r0*np.pi/180.)
    x01=np.linspace(3.,9.5,1300)
    x02=np.linspace(3.,9.5,1300)
    y01=np.interp(x01,x1,diff1)
    y02=np.interp(x02,x2,diff2)
#
#  note the new step size is Delta-x = 6500/1300 pc = 5 pc
#  
    plt.plot(x1,diff1,'r+',label="QI")
    plt.plot(x2,diff2,'b+',label="QIV")
    plt.plot(x01,y01,'r-',label="QI interp")
    plt.plot(x02,y02,'b-',label="QIV interp")
#    plt.plot(x01,y01,'ro',label="QI interp")
#    plt.plot(x02,y02,'bo',label="QIV interp")
    plt.xlabel(r"$s \,(kpc)$")
    plt.ylabel(r"$\Delta \Theta \, (km\, s^{-1})$")
#    plt.legend(loc=2)
    boxx=[3.,3.,9.5,9.5,3.]
    boxy=[-14.5,+14.5,+14.5,-14.5,-14.5]
    plt.plot(boxx,boxy,'k-')
    plt.show()
#
    z=correlate(y01,y02)/(len(y01)*np.std(y01)*np.std(y02))
    z1=correlate(y01,y01)/(len(y01)*np.var(y01))
    z2=correlate(y02,y02)/(len(y01)*np.var(y02))
# just for a test:
##    z=jdccf_0pad(y02,y01)/(np.std(y01)*np.std(y02))
##    z1=jdccf_0pad(y01,y01)/(np.var(y01))
##    z2=jdccf_0pad(y02,y02)/(np.var(y02))
# test done
#    print len(z)
    zx1=np.linspace(1,len(z),len(z))
    zx2=(zx1-(len(z)/2.))*.005
    plt.figure(2)
    plt.plot(zx2,z,'g-',label="ccf")
    plt.plot(zx2,z1,'r-',label="acf Q1")
    plt.plot(zx2,z2,'b-',label="acf Q4")
    plt.legend(loc=1)
##    plt.xlabel(r"$R \,(kpc)$")
##    plt.ylabel(r"$\Delta\Theta (km\, s^{-1})$")
    plt.xlabel(r"$\Delta x \, \, \, (kpc)$")
    plt.ylabel("Normalised Correlation Coefficient")
#    plt.xlim(-3.5,3.5)
#
#  the y axis is normalized to 1 for 100% correlation
#  the x axis is in lag steps of 7.2 pc
#    


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
    return z1,z2,z





def jd_rot_fit2(velsfile1="Q1v.dat",velsfile2="Q4v.dat"):
#
#  this version uses r as the spatial variable, not x (LSCP length)
#
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.signal import medfilt,correlate

    rc('mathtext', default='regular')

    plt.ion()
    plt.clf()
    if velsfile1 is None:
        raise ArgumentError("Please give input file")
    r0 = 8.5
    theta0 = 220.0
    l1,vel1 = get_xy_data(velsfile1)
    l2,vel2 = get_xy_data(velsfile2)
    v1_temp = new_vLSR(l1,vel1)    # Apply the VLSR correction
    v2_temp = new_vLSR(l2,vel2)
    v1 = medfilt(v1_temp,kernel_size=5)
    v2 = medfilt(v2_temp,kernel_size=5)   # Do the median filter that was done in McG07
    ndat1=len(v1)
    ndat2 = len(v2)
    sinl1 = abs(np.sin(l1*np.pi/180))
    wt1= np.array(sinl1)
    sinl2 = abs(np.sin(l2*np.pi/180))
    wt2=np.array(sinl2)
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
#
#   diff1 and diff2 are the residual velocities after subtracting the best fits
#
#  now go back to longitude for equal step in x = linear dist along magic circle
#  we'll just take the range x=3 to x=9.5 kpc (scaled by r0) defined by x01,x02
#
###    x0=abs(l1*r0*np.pi/180.)
    x0=sinl1*r0
#  the 1st quadrant data are in the opposite order now, flip both arrays
    diff0 = diff1
    x1=x0[::-1]
    diff1=diff0[::-1]
###    x2=abs((360.-l2)*r0*np.pi/180.)
    x2=sinl2*r0
###    x01=np.linspace(3.,9.5,1300)
###    x02=np.linspace(3.,9.5,1300)
    x01=np.linspace(3.,7.65,930)
    x02=np.linspace(3.,7.65,930)
    y01=np.interp(x01,x1,diff1)
    y02=np.interp(x02,x2,diff2)
#
#  note the new step size is Delta-x = 6500/1300 pc = 5 pc
#                            Delta-x = 4800/960 pc = 5 pc
#  
    plt.plot(x1,diff1,'r+',label="QI")
    plt.plot(x2,diff2,'b+',label="QIV")
    plt.plot(x01,y01,'r-',label="QI interp")
    plt.plot(x02,y02,'b-',label="QIV interp")
#    plt.plot(x01,y01,'ro',label="QI interp")
#    plt.plot(x02,y02,'bo',label="QIV interp")
    plt.xlabel(r"$r_G \,(kpc)$")
###    plt.xlabel(r"$x \,(kpc)$")
    plt.ylabel(r"$\Delta \Theta \, (km\, s^{-1})$")
#    plt.legend(loc=2)
###    boxx=[3.,3.,9.5,9.5,3.]
    boxx=[3.,3.,7.65,7.65,3.]
    boxy=[-14.5,+14.5,+14.5,-14.5,-14.5]
    plt.plot(boxx,boxy,'k-')
    plt.show()
#
    z=correlate(y01,y02)/(len(y01)*np.std(y01)*np.std(y02))
    z1=correlate(y01,y01)/(len(y01)*np.var(y01))
    z2=correlate(y02,y02)/(len(y01)*np.var(y02))
# just for a test:
#    z=jdccf_0pad(y01,y02)/(np.std(y01)*np.std(y02))
#    z1=jdccf_0pad(y01,y01)/(np.var(y01))
#    z2=jdccf_0pad(y02,y02)/(np.var(y02))
# test done
#    print len(z)
    zx1=np.linspace(1,len(z),len(z))
    zx2=(zx1-(len(z)/2.))*.005
    plt.figure(2)
    plt.plot(zx2,z,'g-',label="ccf")
    plt.plot(zx2,z1,'r-',label="acf Q1")
    plt.plot(zx2,z2,'b-',label="acf Q4")
    plt.legend(loc=1)
##    plt.xlabel(r"$R \,(kpc)$")
##    plt.ylabel(r"$\Delta\Theta (km\, s^{-1})$")
###    plt.xlabel(r"$\Delta x \, \, \, (kpc)$")
    plt.xlabel(r"$\Delta r_G \, (kpc)$")
    plt.ylabel("Normalised Correlation Coefficient")
    plt.minorticks_on()
    vlines(0.,-0.5,1.5,linestyles='dotted',linewidth=1.2)
    #hlines(0.,-3.,3.,linestyles='dotted',linewidth=1.2)
    ylim(-0.5,1.2)
    xlim(-3,3)
##    plt.xlim(-3.5,3.5)
#
#  the y axis is normalized to 1 for 100% correlation
#  the x axis is in lag steps of 7.2 pc
#    


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
    return z1,z2,z
