import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
from scipy.signal import find_peaks
import scipy as sp 

with open('/Users/kartikaarya/Downloads/Internship/CDL_data_41397.pkl',mode='rb') as f:
    data = pkl.load(f)

lam, time, spec = data

fig,ax = plt.subplots(4)

ax[0].set_xlabel('lambda (A)')
ax[0].set_ylabel('Signal (photons/s/steradian/m^2)')
start_time=int(float(input('Start_time(1 - 10s):'))*100)+1
time_slices=int(input(f'Number of time slices that need to be plotted(0 - {(899-start_time)}): '))
ax[0].plot(lam[0:512],spec[0,0:512,start_time:(start_time+time_slices)])
new_lam=lam[0:512]
new_spec=np.empty((time_slices,512),float)

i=0
while i<time_slices:
    new_spec[i]=spec[0,0:512,(start_time+i)]
    i=i+1
        
def gauss_fit3(x, height1, height2, height3, std1, std2, std3, mean1, mean2, mean3): 
    return  height1*np.exp(-0.5*((x-mean1)/std1)**2) + height2*np.exp(-0.5*((x-mean2)/std2)**2) + height3*np.exp(-0.5*((x-mean3)/std3)**2)
x=0
y=0 
avg=0

width1=np.ndarray([])
width1=np.delete(width1,0)
width2=np.ndarray([])
width2=np.delete(width2,0)
width3=np.ndarray([])
width3=np.delete(width3,0)

width1_unc=np.ndarray([])
width1_unc=np.delete(width1_unc,0)
width2_unc=np.ndarray([])
width2_unc=np.delete(width2_unc,0)
width3_unc=np.ndarray([])
width3_unc=np.delete(width3_unc,0)

shift1=np.ndarray([])
shift1=np.delete(shift1,0)
shift1_unc=np.ndarray([])
shift1_unc=np.delete(shift1_unc,0)

while x<time_slices:
    params, pcov= curve_fit(gauss_fit3, new_lam, new_spec[x], p0=(4000000.0, 500000.0, 2000000, 0.04 , 0.04 , 0.04 , 282.1 , 282.40 , 282.55), bounds=((-np.inf,400000.0,-np.inf,-np.inf,0.032,-np.inf,-np.inf,-np.inf,-np.inf,),(np.inf,600000,np.inf,np.inf,0.048,np.inf,np.inf,np.inf,np.inf,)))
    std_all=np.sqrt(np.abs(np.diag(pcov))) # stores the standard deviation in each parameter from the original data
 
    width1= np.insert(width1, x,(params[3]))
    width1_unc= np.insert(width1_unc, x,std_all[3])
    shift1=np.insert(shift1,x, (params[6]-280.448 ))
    shift1_unc= np.insert(shift1_unc, x,std_all[7])
    
    width2= np.insert(width2, x,np.abs((params[4])))
    width2_unc= np.insert(width2_unc, x,(2*std_all[4]))
    
    width3= np.insert(width3, x,np.abs((params[5])))
    width3_unc= np.insert(width3_unc, x,std_all[5])
    
    x+=1
    
    ax[0].plot(lam, gauss_fit3(lam,*params))


ax[1].set_xlabel('time (s)')
ax[1].set_ylabel('width (A)') # plotting standard deviation of each gaussian and how each varies with time. 

ax[1].errorbar(time[start_time:(start_time+time_slices)], width1[0:(time_slices)],yerr=width1_unc[0:(time_slices)], fmt='r.', label='width_peak1')
ax[1].errorbar(time[start_time:(start_time+time_slices)], width2[0:(time_slices)],yerr=width2_unc[0:(time_slices)], fmt='g.', label='width_peak2')
ax[1].errorbar(time[start_time:(start_time+time_slices)], width3[0:(time_slices)],yerr=width3_unc[0:(time_slices)], fmt='b.', label='width_peak3')
leg = ax[1].legend();
ax[1].legend(frameon=True, loc='upper center', ncol=3)
# plotting the error- sqrt of the diagonal of the pcov matrix. for standard deviation of each gaussian fit. 

#converting standard deviation of gaussian fit and the rest wavelength from Å to m

m1=0
lam0_peak1=(280.448)

mc2=1.79*(10**-27)*(3*(10**8))**2
c=3*(10**8)

a=0

temp_peak1=np.ndarray([])
temp_peak1=np.delete(temp_peak1,0)
temp_unc1=np.ndarray([])
temp_unc1=np.delete(temp_unc1,0)
unc1=0
temp1=0
e=1.6*10**-19
while y<time_slices:
    temp1=((((width1[y]**2)/(lam0_peak1**2))*mc2 )/e)
    temp_peak1=np.insert(temp_peak1,y,temp1)
    unc1=(2*mc2*(width1_unc[y]**2))/(lam0_peak1**2)
    temp_unc1=np.insert(temp_peak1,y,temp1)
    y+=1
ax[2].errorbar(time[start_time:(start_time+time_slices)], temp_peak1[0:(time_slices)],yerr=temp_unc1[0:(time_slices)], fmt='m.', label='Temperature calculated from peak 1')
ax[2].set_xlabel('time (s)')
ax[2].set_ylabel('Energy (eV)') # plotting standard deviation of each gaussian and how each varies with time.
leg1 = ax[2].legend();
ax[2].legend(frameon=True, loc='upper center', ncol=1,fancybox=True, framealpha=1)

speed_peak1=np.ndarray([])
speed_peak1=np.delete(speed_peak1,0)

speed=0
z=0
speed_unc1=np.ndarray([])
speed_unc1=np.delete(speed_unc1,0)

sunc1=0
while z<time_slices:
    speed=((shift1[z])/(lam0_peak1))*c 
    speed_peak1=np.insert(speed_peak1,z,speed)
    sunc1=((shift1_unc[z])/(lam0_peak1))*c 
    speed_unc1=np.insert(speed_peak1,z,sunc1)
    z+=1
    
ax[3].errorbar(time[start_time:(start_time+time_slices)],speed_peak1[0:(time_slices)],yerr=speed_unc1[0:(time_slices)], fmt='c.', label='speed calculated from shift in peak 1')
ax[3].set_xlabel('time (s)')
ax[3].set_ylabel('speed (m/s)') 
leg2 = ax[3].legend();
ax[3].legend(frameon=True, loc='upper center', ncol=1,fancybox=True, framealpha=1)
print (speed_unc1, temp_unc1 )


print(f'''lam {lam.shape}
time {time.shape}
spec {spec.shape}
new_lam {new_lam.shape}
new_spec {new_spec.shape}
width1 {width1.shape}
width2 {width2.shape}
width3 {width3.shape}
width1_unc {width1_unc.shape}
width2_unc {width2_unc.shape}
width3_unc {width3_unc.shape}
shift1 {shift1.shape}
shift1_unc {shift1_unc.shape}
std_all {std_all.shape}
temp_peak1 {temp_peak1.shape}
temp_unc1 {temp_unc1.shape}
unc1 {unc1.shape}
temp1 {temp1.shape}
speed_peak1 {speed_peak1.shape}
speed {speed.shape}
speed_unc1 {speed_unc1.shape}''')
plt.show() 