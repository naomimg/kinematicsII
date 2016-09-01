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
    Utility function to retrieve x,y vectors form an ASCII data file.
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
    Uo = 14.60                # km/s
    Vo =  9.6
    Wo =  9.3
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

   
### Bits for plotting Clemens data
clem_data=ascii.read("clemens85.txt", comment="#")  
#    dat = read_data("rotvels.dat",(0,1))
c_lon = data['col1']
c_vel = data['col2']


dat = ascii.read("fitpars.dat")
l1 = dat['col1']
l2 = dat['col2']
vo = dat['col4']



# Bits for plotting the two fitpars
twocomp_data=ascii.read("outspec_2comp.dat", comment="#")  
vel2 = twocomp_data['col1']
Tb2 = twocomp_data['col2']
Tbfit2 = twocomp_data['col3']
semilogy(vel2,Tb2,'ks',mfc='none')
semilogy(vel2,Tbfit2,'k--',label='2 component')

threecomp_data=ascii.read("outspec_3comp.dat", comment="#")  
vel3 = threecomp_data['col1']
Tb3 = threecomp_data['col2']
Tbfit3 = threecomp_data['col3']
semilogy(vel3,Tb3,'ks',mfc='none')
semilogy(vel3,Tbfit3,'k-',label='3 component')

ylabel(r'$T_{b}\, (K)$')
xlabel(r'$V_{r}\, (km\, s^{-1})$')
xlim(-20.,50.)
ylim(0.01,150.)
legend()



# Bits for resampling the rotation curve and doing ACF analysis

velsfile1="firstquad_vels.dat"
velsfile2="../fourthquad_vels.dat"
r0 = 8.5
theta0 = 220.0
#    r0=8.34    # Reid et al (2014) values
#    theta0=240.0
l1,vel1 = get_xy_data(velsfile1)
l2,vel2 = get_xy_data(velsfile2)
#v1 = vel1
#v2 = medfilt(vel2,kernel_size=5)
v1 = new_vLSR(l1,vel1)    # Apply the VLSR correction
v2_temp = new_vLSR(l2,vel2)
v2 = medfilt(v2_temp,kernel_size=5)

ndat1=len(v1)
ndat2 = len(v2)
wt1= array(sinl1)
wt2=array(sinl2)

sinl1 = abs(sin(l1*pi/180))
sinl2 = abs(sin(l2*pi/180))
rot1 = abs(v1) + theta0*sinl1
rg1 = r0*sinl1
rot2 = abs(v2) + theta0*sinl2
rg2 = r0*sinl2

# Create arrays with both datasets    rot1 = array(abs(v1) + theta0*sinl1)
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

# Create arrays of the differences from fit
fit1 = polyval(pars,rg1/r0)
fit2 = polyval(pars,rg2/r0)
diff1 = rot1 - fit1*theta0
diff2 = rot2 - fit2*theta0



# Get the data
velsfile1="firstquad_vels.dat"
velsfile2="../fourthquad_vels.dat"
l1,vel1_temp = get_xy_data(velsfile1)
l2,vel2_temp = get_xy_data(velsfile2)
vel1 = medfilt(vel1_temp,kernel_size=15)
vel2 = medfilt(vel2_temp,kernel_size=15)


# Resample the data based on a common l
f1=interp1d(l1,vel1)
l1_r = np.linspace(22.0,65.0,550)
vel1_r=f1(l1_r)
f2=interp1d(abs(l2-360.),-vel2)
vel2_r=f2(l1_r)
diff_vel = vel1_r - vel2_r
plot(l1_r,diff_vel,'k-')
xlabel(r"$\left|Galactic\, Longitude\, (deg)\right|$")
ylabel(r"$\Delta V_t\, (km\,s^{-1})$")
m=mean(diff_vel[420:])
print "Mean V_diff old LSR: %.2f" % m
#savefig("vdiff_vs_l.pdf")

# Do the same thing with newVLSR values
v1 = new_vLSR(l1,vel1)    # Apply the VLSR correction
v2 = new_vLSR(l2,vel2)
#v2 = medfilt(v1_temp,kernel_size=15)
#v2 = medfilt(v2_temp,kernel_size=15)
f1=interp1d(l1,v1)
l1_r = np.linspace(22.0,65.0,550)
vel1_r=f1(l1_r)
f2=interp1d(abs(l2-360.),-v2)
vel2_r=f2(l1_r)
diff_vel = vel1_r - vel2_r
plot(l1_r,diff_vel,'k--')
m=mean(diff_vel[420:])
print "Mean V_diff new LSR: %.2f" % m
xlabel(r"$\left| Galactic\, Longitude\, (deg)\right|$")
ylabel(r"$\Delta V_t\, (km\,s^{-1})$")

# Overplot Blitz \& Spergel (1991)
a = [7.,10.,14.]
l_mod = arange(0.,180.)
for val in a:
    deltav = 2.* val * np.cos(l_mod*math.pi/180.)
    plot(l_mod,deltav,'k:')
xlim(20,65.)
ylim(-20.,25.)
minorticks_on()
pylab.text(40,21.4,'$\Pi_0=14$')
pylab.text(25,18.3,'$\Pi_0=10$')
pylab.text(23,13,'$\Pi_0=7$')
savefig("vdiff_vs_l.pdf")


# Resample the data based on rg
f1=interp1d(rg1,diff1)
rg1_r = np.linspace(3.0,max(rg1),550)
diff1_r=f1(rg1_r)
# Resample the RG data
f2=interp1d(rg2,diff2)
diff2_r=f2(rg1_r)
diff_rot = diff1_r - diff2_r
plot(rg1_r,diff_rot,'k-')
xlabel(r"$R\, (kpc)$")
ylabel(r"$\Delta \Theta\, (km\,s^{-1})$")
#savefig("vdiff_vs_r.pdf")

# Power spectrum
data =diff1_r
ps = np.abs(np.fft.fft(data))**2
r_step = abs(rg1_r[1]-rg1_r[0])
freqs = np.fft.fftfreq(data.size, r_step)
idx = np.argsort(freqs)
plt.semilogy(freqs[idx], ps[idx])
xlabel(r'$k\, (kpc^{-1})$')

# Auto-correlation function
x=diff1_r
x = x - x.mean()
autocorr = np.correlate(x, x, mode='full')
autocorr = autocorr[x.size:]
autocorr /=autocorr.max()
plot(autocorr,label='ACF QI')
x=diff2_r
x = x - x.mean()
autocorr = np.correlate(x, x, mode='full')
autocorr = autocorr[x.size:]
autocorr /=autocorr.max()
plot(autocorr,label='ACF QIV')

# Cross correlation
ccf=np.correlate(diff1_r,diff2_r,'full')
ccf = ccf[x.size:]
ccf /=ccf.max()
plot(ccf,label='CCF')
xlabel('Lags')
ylabel('Correlation')

#
f=rot1_r
acf=np.correlate(f,f)

FFTacf=fft(acf);
PSacf=abs(FFTacf);


# Also do the resampling in l for a difference in l curve
f1=interp1d(l1,v1)
l_r = np.linspace(20.5,66.,550)
y2= abs(l2-360.)
f2=interp1d(y2,abs(v2))
v1_r = f1(l_r)
v2_r = f2(l_r)
v_diff = v1_r - v_2
plot(l_r,v_diff,'k-')
xlabel(r'$Galactic\, Longitude\, (deg)$')
ylabel(r'$\Delta V\, (km\, s^{-1})$')
savefig("vdiff_vs_l.pdf")





plot(rg1,rot1,'r+',label="QI")
plot(rg2,rot2,'b+',label="QIV")
minorticks_on()
legend()
xlabel(r"$R\, (kpc)$")
ylabel(r"$\Theta\, (km\,s^{-1})$")
print 'Number of 1st quad points: ',len(rot1) 
print 'Number of 1st quad points: ',len(rot2)
xlim(3.0,8.0)
ylim(190.,260.)



# Universal rotation curve
r0=8.34
a1=241 
a2 = 0.9 
a3= 1.46 
#beta = 1.10 * G * M / (a1**2 * ropt)

for a2 in arange(0.6,1.2,0.1):
	ropt = a2* r0
	x = urg/ropt
	Vdsq = a1**2 * a2 * (1.97*x**1.22)/(x**2+0.78**2)**1.43
	Vhsq = a1**2 * (1-a2)*(1+a3**2)*x**2/(x**2 + a3**2)

	ll = a3
	V = a1 *(((0.72 + 0.44*log10(ll))*(1.97*x**1.22)/(x**2+0.78**2)**1.43) + 1.6*exp(-0.4*ll)*x**2/(x**2 + 1.5**2 * ll**0.4))

	#V = (Vdsq + Vhsq)**0.5
	lstr = '%.1f' %a2
	plot(urg,V,label=lstr)
legend()

ll = a3
V = a1 *(((0.72 + 0.44*log10(ll))*(1.97*x**1.22)/(x**2+0.78**2)**1.43) + 1.6*exp(-0.4*ll)*x**2/(x**2 + 1.5**2 * ll**0.4))


# Get the fraction of pixels above the masking limit (40.0 K)
import glob
fracs={}
mosfiles=glob.glob("*MOS*")
for filename in mosfiles:
    lats=[]
    hdulist=fits.open(filename)
    im=hdulist[0].data
    dimensions = im.shape
    hdr=hdulist[0].header
    idimensions=0
    if dimensions[0]==1:
        idimensions=1        ## this means that the header most likely has a degenerate polarization axis
    nypix=dimensions[1+idimensions]
    nxpix=dimensions[2+idimensions]

    lat_ref_pix = hdr['crpix2']
    lat_del_pix = hdr['cdelt2']
    lat_ref_value = hdr['crval2']   
    lats=np.array(hdr['crval2']+hdr['cdelt2']*(np.arange(nypix) + 1 - hdr['crpix2']))
    count=0
    for i in range(dimensions[2]):
        for j in range(dimensions[3]):
            if lats[i]>=-0.5 and lats[i]<=0.5 and im[0,0,i,j]>=30.0:
                count+=1
    print count
    fracs[filename]=float(count)/(float(dimensions[2])*float(dimensions[3]))*100.

# Get the average value
av=0
for key in fracs.keys():
    print fracs[key]
    av=av+fracs[key]

print av/len(fracs)






