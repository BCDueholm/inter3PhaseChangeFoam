import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('results', delimiter=' ',skiprows=0)

cell_vol = 0.1*0.1*(2/640)
domain_vol = 0.1*0.1*2
rho_l = 1000

################################# GLOBAL ERROR #############################
error_global = (abs((data[len(data[:,1])-1,1] 
            - data[0,1])-(data[len(data[:,2])-1,2]- data[0,2]))/data[0,1])*100
################################ MASS RELATIVE ERROR #########################
error = np.zeros(len(data[:,1])-1)

for i in range(0,len(error)):
    error[i] = (abs((data[i+1,1]-data[i,1])
                    -(data[i+1,2]-data[i,2]))/data[0,1])*100
   
plt.plot(data[:len(error),0],error)
plt.xlabel('Time [s]')
plt.ylabel('Mass relative error [%]')
plt.ylim(0,1e-4)
plt.xlim(0,0.1)
plt.grid(axis='both')
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.savefig('MassRelativeError.png')

plt.close()

############################ SUM OF ALPHAS ##################################

avg_sum_alpha = (((data[:,1]/rho_l)+data[:,2]+data[:,3])/cell_vol)/640
plt.plot(data[:,0],avg_sum_alpha)
plt.ylim(0,2)
plt.xlim(0,0.1)
plt.xlabel(r'Time [s]')
plt.ylabel(r'Average sum of void fractions [-]')
plt.grid(axis='both')
plt.savefig('AverageSumOfAlphas.png')
plt.close()

#################### Volume average development of alphas ####################

alpha_l = (data[:,1]/rho_l)/domain_vol
alpha_v = data[:,2]/domain_vol
alpha_nc = data[:,3]/domain_vol

plt.plot(data[:,0],alpha_l, label ='alpha_l', linestyle='-.')
plt.plot(data[:,0],alpha_v, label ='alpha_v', linestyle='--')
plt.plot(data[:,0],alpha_nc, label ='alpha_nc')
plt.ylim(0,1)
plt.xlim(0,0.1)
plt.xlabel(r'Time [s]')
plt.ylabel(r'Volume averaged void fractions [-]')
plt.grid(axis='both')
plt.legend(loc='upper right')
plt.savefig('VolumeAveragedAlphas.png')
plt.close()

############################# Mass alpha_l ###################################

plt.plot(data[:,0],data[:,1])
plt.ylim(9.49,9.5)
plt.xlim(0,0.1)
plt.xlabel(r'Time [s]')
plt.ylabel(r'Mass liquid [kg]')
plt.grid(axis='both')
plt.savefig('MassOfLiquid.png')
plt.close()

# ############################# Mass alpha_v ###################################

plt.plot(data[:,0],data[:,2])
plt.ylim(0,1e-2)
plt.xlim(0,0.1)
plt.xlabel(r'Time [s]')
plt.ylabel(r'Mass vapour [kg]')
plt.grid(axis='both')
plt.savefig('MassOfVapour.png')
plt.close()

# ############################# Mass alpha_nc ###################################

plt.plot(data[:,0],data[:,3])
plt.ylim(0,1e-2)
plt.xlim(0,0.1)
plt.xlabel(r'Time [s]')
plt.ylabel(r'Mass air [kg]')
plt.grid(axis='both')
plt.savefig('MassOfAir.png')
plt.close()

############################# Surface height lv ###############################

plt.plot(data[:,0],data[:,4])
plt.ylim(0.92,0.96)
plt.xlim(0,0.1)
plt.xlabel('Time [s]')
plt.ylabel('Liquid/vapour surface height [m]')
plt.grid(axis='both')
plt.savefig('SurfaceHeightLV.png')
plt.close()

# ############################# Surface height vnc #############################

plt.plot(data[:,0],data[:,5])
plt.ylim(1,1.8)
plt.xlim(0,0.1)
plt.xlabel('Time [s]')
plt.ylabel('Vapour/air surface height [m]')
plt.grid(axis='both')
plt.savefig('SurfaceHeightVNC.png')
plt.close()

# ############################# Surface velocity lv #############################

plt.plot(data[:,0],data[:,6])
plt.ylim(0,0.45)
plt.xlim(0,0.1)
plt.xlabel('Time [s]')
plt.ylabel('Liquid/vapour interface velocity [m/s]')
plt.grid(axis='both')
plt.savefig('SurfaceVelocityLV.png')
plt.close()

# ############################# Surface velocity vnc #############################

plt.plot(data[:,0],data[:,7])
plt.ylim(0,13)
plt.xlim(0,0.1)
plt.xlabel('Time [s]')
plt.ylabel('Vapour/air interface velocity [m/s]')
plt.grid(axis='both')
plt.savefig('SurfaceVelocityVNC.png')

plt.close()