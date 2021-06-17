# Importing the ChoiceMC class, this structure should be improved
import sys
try:
    from ChoiceMC import ChoiceMC
except ModuleNotFoundError:
    sys.path.append('..')
    from ChoiceMC import ChoiceMC
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import time
import os

# Making a folder to store the current test
time_str = "NMM_Sweep-"+str(time.gmtime().tm_year)+'-'+str(time.gmtime().tm_mon)+'-'+str(time.gmtime().tm_mday)
path = os.path.join(os.getcwd(), time_str)
try:
    os.mkdir(path)
except FileExistsError:
    pass
os.chdir(path)

# Making a folder to store the individual results from each data point
data_path = os.path.join(path, "Individual_Results")
try:
    os.mkdir(data_path)
except FileExistsError:
    pass

# Setting up the variables to sweep over
g_sweep = [0.1, 1.0, 2.0, 4.0]
mMax_sweep = range(1,26)
P_sweep = range(1,42,2)

g_sweep_dict = {}
for ig, g in enumerate(g_sweep):
    energy = np.zeros((len(P_sweep),2), float)
    energy_out = open(os.path.join(data_path,"EnergyPSweep_g" + str(g) + '.dat'),'w')
    for iP, P in enumerate(P_sweep):
        PIMC = ChoiceMC(m_max=25, P=P, g=g, MC_steps=10000, N=2, PIGS=True, Nskip=100, Nequilibrate=100)
        PIMC.createRhoNMM()
        energy[iP,:] = [P, PIMC.E0_nmm]
        energy_out.write(str(P) + ' ' + str(PIMC.E0_nmm) + '\n')
    
    g_sweep_dict.update({g: energy})
    
    # Plotting
    E_fig, E_ax = plt.subplots(1, 1, figsize=(8,5))
    E_ax.plot(energy[:,0], energy[:,1])
    E_ax.set_xlabel('P')
    E_ax.set_ylabel(r'$E_0$')
    E_ax.annotate('N = 2; g = ' + str(g) + r'$;\ M_{Max}$' + ' = 25', xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
    E_ax.minorticks_on()
    E_fig.tight_layout()
    E_fig.savefig(os.path.join(data_path,"EnergyPSweep_g" + str(g) + ".png"))
    
    energy_out.close()
    plt.close('all')

# Plotting the energy versus g for varied mMax
E_fig = plt.figure(figsize=(8,3*((len(g_sweep)+1)//2)))
for i, g in enumerate(g_sweep_dict):
    ax = E_fig.add_subplot((len(g_sweep)+1)//2, 2, i+1)
    ax.plot(g_sweep_dict[g][:,0], g_sweep_dict[g][:,1])
    ax.annotate('N = 2; g = ' + str(g) + r'$;\ M_{Max}$' + ' = 25', xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
    ax.set_xlabel('P')
    ax.set_ylabel(r'$E_0$')
    ax.minorticks_on()
E_fig.tight_layout()
E_fig.savefig("Energy_PSweep.png")
plt.close('all')

g_sweep_dict = {}
for ig, g in enumerate(g_sweep):
    energy = np.zeros((len(mMax_sweep),2), float)
    energy_out = open(os.path.join(data_path,"EnergymMaxSweep_g" + str(g) + '.dat'),'w')
    for imMax, mMax in enumerate(mMax_sweep):
        PIMC = ChoiceMC(m_max=mMax, P=41, g=g, MC_steps=10000, N=2, PIGS=True, Nskip=100, Nequilibrate=100)
        PIMC.createRhoNMM()
        energy[imMax,:] = [mMax, PIMC.E0_nmm]
        energy_out.write(str(mMax) + ' ' + str(PIMC.E0_nmm) + '\n')
    
    g_sweep_dict.update({g: energy})
    
    # Plotting
    E_fig, E_ax = plt.subplots(1, 1, figsize=(8,5))
    E_ax.plot(energy[:,0], energy[:,1])
    E_ax.set_xlabel(r'$M_{Max}$')
    E_ax.set_ylabel(r'$E_0$')
    E_ax.annotate('N = 2; g = ' + str(g) + '; P = 41', xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
    E_ax.minorticks_on()
    E_fig.tight_layout()
    E_fig.savefig(os.path.join(data_path,"EnergymMaxSweep_g" + str(g) + ".png"))
    
    energy_out.close()
    plt.close('all')

# Plotting the energy versus g for varied mMax
E_fig = plt.figure(figsize=(8,3*((len(g_sweep)+1)//2)))
for i, g in enumerate(g_sweep_dict):
    ax = E_fig.add_subplot((len(g_sweep)+1)//2, 2, i+1)
    ax.plot(g_sweep_dict[g][:,0], g_sweep_dict[g][:,1])
    ax.annotate('N = 2; g = ' + str(g) + '; P = 41', xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
    ax.set_xlabel(r'$M_{Max}$')
    ax.set_ylabel(r'$E_0$')
    ax.minorticks_on()
E_fig.tight_layout()
E_fig.savefig("Energy_mMaxSweep.png")
plt.close('all')

    
