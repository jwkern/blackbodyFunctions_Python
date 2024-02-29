#import matplotlib
#matplotlib.use('TkAgg')


import matplotlib.pyplot as plt
import math
import numpy as np
import pylab

##################################################################################
######### Stellar Arrays & Constants #############################################
##################################################################################

T = 3000                                #stellar temperature (K)
wl_steps = 150000                       #wavelength steps
wl = np.linspace(0, 0.01, wl_steps)     #wavelength array
h = 6.626e-34                           #plancks constant (m^2*kg*s^-1)
c = 3e8                                 #speed of light (m*s^-1)
k = 1.38e-23                            #boltzmann const. (m^2*kg*s^-2*K^-1)
e = 2.71                                #exponential
Bs = np.zeros(wl_steps)                 #blackbody function array
Bsnorm = np.zeros(wl_steps)             #normalized blackbody value
Bg1 = np.zeros(wl_steps)                #blackbody function array
Bg2 = np.zeros(wl_steps)                #blackbody function array
Bg1norm = np.zeros(wl_steps)            #normalized blackbody value
Bg2norm = np.zeros(wl_steps)            #normalized blackbody value
Btotal = np.zeros(wl_steps)             #blackbody function array
Btotalnorm = np.zeros(wl_steps)         #normalized blackbody value
Bg1_total = np.zeros(wl_steps)          #total blackbody from dust 1
Bg2_total = np.zeros(wl_steps)          #total blackbody from dust 2



###################################################################################
######### Stellar blackbody function 1 ############################################
###################################################################################

for i in range (1, wl_steps):
    Bs[i] = ((2.0*h*c**2)/(wl[i]**5))/(e**((h*c)/(wl[i]*k*T)) - 1.0)

Bsmax = Bs.max()

for i in range (1, wl_steps):
    Bsnorm[i] = Bs[i] / Bsmax


##################################################################################
######### Dust grain blackbody function 1 ########################################
##################################################################################
Bg1 = np.zeros(wl_steps)                #blackbody function array
Bg1norm = np.zeros(wl_steps)            #normalized blackbody value
Bg1_total = np.zeros(wl_steps)          #total blackbody from dust 1

sun_radius = 0.00465046324              #sun radius in AU
d1 = 0.2                                 #distance in AU
frac_rs_d1 = sun_radius/d1                #ratio of stars radius to distance from star
Tg1=(T/2**(0.5))*(frac_rs_d1)**(0.5)     #temperature of dust grain
particle_number1=1e3                     #number of particles
print('Disk 1 distance = ' + str(d1))
print('Disk 1 temperature = ' + str(Tg1))
print('Disk 1 particle number = ' + str(particle_number1))

for i in range (1, wl_steps):
    Bg1[i] = ((2.0*h*c**2)/(wl[i]**5))/(e**((h*c)/(wl[i]*k*Tg1)) - 1.0)
    Bg1_total[i] = particle_number1*Bg1[i]

Bg1max = Bg1.max()

for i in range (1, wl_steps):
    Bg1norm[i] = Bg1[i] / Bg1max


##################################################################################
######### Dust grain blackbody function 2 ########################################
##################################################################################
Bg2 = np.zeros(wl_steps)                #blackbody function array
Bg2norm = np.zeros(wl_steps)            #normalized blackbody value
Bg2_total = np.zeros(wl_steps)          #total blackbody from dust 2

d2 = 0.1                                 #distance in AU
frac_rs_d2 = sun_radius/d2                #ratio of stars radius to distance from star
Tg2=(T/2**(0.5))*(frac_rs_d2)**(0.5)     #temperature of dust grain
particle_number2=2e3                     #number of particles
print ('Disk 2 distance = ' + str(d2))
print ('Disk 2 temperature = ' + str(Tg2))
print ('Disk 2 particle number = ' + str(particle_number2))

for i in range (1, wl_steps):
    Bg2[i] = ((2.0*h*c**2)/(wl[i]**5))/(e**((h*c)/(wl[i]*k*Tg2)) - 1.0)
    Bg2_total[i] = particle_number2*Bg2[i]

Bg2max = Bg2.max()

for i in range (1, wl_steps):
    Bg2norm[i] = Bg2[i] / Bg2max


##################################################################################
######### Dust grain blackbody function 2 ########################################
##################################################################################
Bg3 = np.zeros(wl_steps)                #blackbody function array
Bg3norm = np.zeros(wl_steps)            #normalized blackbody value
Bg3_total = np.zeros(wl_steps)          #total blackbody from dust 2

d3 = 0.05                                #distance in AU
frac_rs_d3 = sun_radius/d3                #ratio of stars radius to distance from star
Tg3=(T/2**(0.5))*(frac_rs_d3)**(0.5)     #temperature of dust grain
particle_number3=3e3                     #number of particles
print('Disk 3 distance = ' + str(d3))
print('Disk 3 temperature = ' + str(Tg3))
print('Disk 3 particle number = ' + str(particle_number3))



for i in range (1, wl_steps):
    Bg3[i] = ((2.0*h*c**2)/(wl[i]**5))/(e**((h*c)/(wl[i]*k*Tg3)) - 1.0)
    Bg3_total[i] = particle_number3*Bg3[i]

Bg3max = Bg3.max()

for i in range (1, wl_steps):
    Bg3norm[i] = Bg3[i] / Bg3max




###################################################################################
######### Combined blackbody function #############################################
###################################################################################
wl_microns = np.zeros(wl_steps)

for i in range (1, wl_steps):
    Btotal[i] = Bs[i] + Bg1_total[i] + Bg2_total[i] + Bg3_total[i] 

Btotalmax = Btotal.max()

for i in range (1, wl_steps):
    Btotalnorm[i] = Btotal[i] / Btotalmax
    wl_microns[i] = wl[i]*1e6


###################################################################################
##########   Plot Details #########################################################
###################################################################################

plt.plot (wl_microns, Bs, '-.', label = 'Star')
plt.plot (wl_microns, Bg1_total, '--', label = 'Disk 1')
plt.plot (wl_microns, Bg2_total, '--', label = 'Disk 2')
plt.plot (wl_microns, Bg3_total, '--', label = 'Disk 3')
plt.plot (wl_microns, Btotal, '-', label = 'Combined')
plt.legend (loc='upper right')
plt.xlabel('$\lambda$ ($\mu$m)')
plt.ylabel('B($\lambda$, T)')
plt.title('Blackbody of Exoplanet Disk System')
plt.xlim (-1e-7, 25)
#plt.ylim (0, 1500)
plt.show()
