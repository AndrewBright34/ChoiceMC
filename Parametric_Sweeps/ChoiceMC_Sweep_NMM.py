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

parent_dir = os.getcwd()
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
mMax_sweep = [5, 10, 15, 25]
g_sweep = [0.1, 1.0, 2.0, 4.0]
P_sweep = range(1,42,2)

for mMax in mMax_sweep:
    g_sweep_dict = {}
    for ig, g in enumerate(g_sweep):
        energy = np.zeros((len(P_sweep),2), float)
        energy_out = open(os.path.join(data_path,"EnergyPSweep_mMax" + str(mMax) + "_g" + str(g) + '.dat'),'w')
        for iP, P in enumerate(P_sweep):
            PIMC = ChoiceMC(m_max=mMax, P=P, g=g, MC_steps=10000, N=2, PIGS=True, Nskip=100, Nequilibrate=100, T=2./3.)
            PIMC.runNMM()
            energy[iP,:] = [P, PIMC.E0_NMM]
            energy_out.write(str(P) + ' ' + str(PIMC.E0_NMM) + '\n')
        
        g_sweep_dict.update({g: energy})
        
        # Plotting
        E_fig, E_ax = plt.subplots(1, 1, figsize=(8,5))
        E_ax.plot(energy[:,0], energy[:,1])
        E_ax.set_xlabel('P')
        E_ax.set_ylabel(r'$E_0$')
        E_ax.annotate('N = 2; g = ' + str(g) + r'$;\ M_{Max}$' + ' = ' + str(mMax), xy=(0.5, 0.15),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
        E_ax.minorticks_on()
        E_fig.tight_layout()
        E_fig.savefig(os.path.join(data_path,"EnergyPSweep_mMax" + str(mMax) + "_g" + str(g) + ".png"))
        
        energy_out.close()
        plt.close('all')
    
    # Plotting the energy versus g for varied mMax
    E_fig = plt.figure(figsize=(8,3*((len(g_sweep)+1)//2)))
    for i, g in enumerate(g_sweep_dict):
        ax = E_fig.add_subplot((len(g_sweep)+1)//2, 2, i+1)
        ax.plot(g_sweep_dict[g][:,0], g_sweep_dict[g][:,1])
        ax.annotate('N = 2; g = ' + str(g) + r'$;\ M_{Max}$' + ' = ' + str(mMax), xy=(0.5, 0.15),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
        ax.set_xlabel('P')
        ax.set_ylabel(r'$E_0$')
        ax.minorticks_on()
    E_fig.tight_layout()
    E_fig.savefig("Energy_PSweep_mMax" + str(mMax) + ".png")
    plt.close('all')

# Setting up the variables to sweep over
P_sweep = [3, 13, 21, 51]
g_sweep = [0.1, 1.0, 2.0, 4.0]
mMax_sweep = range(1,26)

for P in P_sweep:
    g_sweep_dict = {}
    for ig, g in enumerate(g_sweep):
        energy = np.zeros((len(mMax_sweep),2), float)
        energy_out = open(os.path.join(data_path,"EnergymMaxSweep_P" + str(P) + "_g" + str(g) + '.dat'),'w')
        for imMax, mMax in enumerate(mMax_sweep):
            PIMC = ChoiceMC(m_max=mMax, P=P, g=g, MC_steps=10000, N=2, PIGS=True, Nskip=100, Nequilibrate=100, T=2./3.)
            PIMC.runNMM()
            energy[imMax,:] = [mMax, PIMC.E0_NMM]
            energy_out.write(str(mMax) + ' ' + str(PIMC.E0_NMM) + '\n')
        
        g_sweep_dict.update({g: energy})
        
        # Plotting
        E_fig, E_ax = plt.subplots(1, 1, figsize=(8,5))
        E_ax.plot(energy[:,0], energy[:,1])
        E_ax.set_xlabel(r'$M_{Max}$')
        E_ax.set_ylabel(r'$E_0$')
        E_ax.annotate('N = 2; g = ' + str(g) + '; P = 9', xy=(0.5, 0.15),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
        E_ax.minorticks_on()
        E_fig.tight_layout()
        E_fig.savefig(os.path.join(data_path,"EnergymMaxSweep_P" + str(P) + "_g" + str(g) + ".png"))
        
        energy_out.close()
        plt.close('all')
    
    # Plotting the energy versus g for varied mMax
    E_fig = plt.figure(figsize=(8,3*((len(g_sweep)+1)//2)))
    for i, g in enumerate(g_sweep_dict):
        ax = E_fig.add_subplot((len(g_sweep)+1)//2, 2, i+1)
        ax.plot(g_sweep_dict[g][:,0], g_sweep_dict[g][:,1])
        ax.annotate('N = 2; g = ' + str(g) + '; P = 9', xy=(0.5, 0.15),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
        ax.set_xlabel(r'$M_{Max}$')
        ax.set_ylabel(r'$E_0$')
        ax.minorticks_on()
    E_fig.tight_layout()
    E_fig.savefig("Energy_mMaxSweep_P" + str(P) + ".png")
    plt.close('all')

# Sweeping over g values - this is what will be used in future tests
g_sweep = np.sort(np.append(np.linspace(0.01, 8, 40),np.array([0.1, 1., 2.])))
energy = np.zeros((len(g_sweep),2), float)
energy_out = open("Energy_NMM.dat", 'w')
orrCorr = np.zeros((len(g_sweep),2), float)
orrCorr_out = open("OrrCorr_NMM.dat", 'w')
for ig, g in enumerate(g_sweep):
    PIMC = ChoiceMC(m_max=10, P=21, g=g, MC_steps=10000, N=2, PIGS=True, Nskip=100, Nequilibrate=100, T=2./3.)
    PIMC.runNMM()
    energy[ig,:] = [g, PIMC.E0_NMM]
    energy_out.write(str(g) + ' ' + str(PIMC.E0_NMM) + '\n')
    orrCorr[ig,:] = [g, PIMC.eiej_NMM]
    orrCorr_out.write(str(g) + ' ' + str(PIMC.eiej_NMM) + '\n')

# Plotting E0 versus g
fig, ax = plt.subplots(1, 1, figsize=(8,5))
ax.plot(energy[:,0], energy[:,1], label='NMM', marker='o', color='k')
ax.set_xlabel('g')
ax.set_ylabel(r'$E_0$')
ax.annotate('N = 2; P = 21;' + r'$\ M_{Max}$' + ' = 10', xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
ax.minorticks_on()
ax.legend()
fig.tight_layout()
fig.savefig("Energy_NMM.png")

# Plotting O vs g
O_fig, O_ax = plt.subplots(1, 1, figsize=(8,5))
O_ax.plot(orrCorr[:,0], orrCorr[:,1], label='NMM', marker='o', color='k')
O_ax.set_xlabel('g')
O_ax.set_ylabel('Orientational Correlation')
ax.annotate('N = 2; P = 21;' + r'$\ M_{Max}$' + ' = 10', xy=(0.5, 0.15),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
O_ax.minorticks_on()
ax.legend()
O_fig.tight_layout()
O_fig.savefig("OrrCorr_NMM.png")

energy_out.close()
orrCorr_out.close()
plt.close('all')

# Performing sweeps over tau and beta
tau_sweep = [0.0015, 0.005, 0.05, 0.1]
g_sweep = [0.1, 1., 2., 4.]
P_sweep = [3, 5, 7, 9, 11, 13, 15, 19, 23, 27, 31, 35]
N = 2
for tau in tau_sweep:
    # Making a folder to store the current test
    time_str = "Beta_Sweep_NMM_tau"+str(tau)+"-"+str(time.gmtime().tm_year)+'-'+str(time.gmtime().tm_mon)+'-'+str(time.gmtime().tm_mday)
    path = os.path.join(parent_dir, time_str)
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

    T_sweep = []
    for P in P_sweep:
        T_sweep.append(1./(tau*P))

    # Creating dictionaries to store the results from the g and N sweep
    g_sweep_dict_E = {}
    # Performing the sweep over mMax
    for ig, g in enumerate(g_sweep):
        print("------------------------------------------------")
        print("Starting g = " + str(g))
        
        # Creating arrays to store the E versus g
        energy = np.zeros((len(T_sweep),2), float)
        energy_out = open(os.path.join(data_path,"Energy_g" + str(g) + '.dat'),'w')
    
        for i, T in enumerate(T_sweep):
            print("------------------------------------------------")
            print("Starting T = " + str(T) + "; P = " + str(P_sweep[i]))
            # Creating a ChoiceMC object for the current iteration
            PIMC = ChoiceMC(m_max=25, P=P_sweep[i], g=g, MC_steps=1000, N=N, PIGS=True, Nskip=100, Nequilibrate=100, T=T)
            print("Beta = " + str(PIMC.beta) + "; Tau = " + str(PIMC.tau))
            # Creating the probability density matrix for each rotor
            PIMC.runNMM()
            
            # Storing and saving the data from the current run
            energy[i,:] = [PIMC.beta, PIMC.E0_NMM]
            energy_out.write(str(PIMC.beta) + ' ' + str(PIMC.E0_NMM) + '\n')
                
            # Closing the remaining open plots
            plt.close('all')
    
        # Plotting
        E_fig, E_ax = plt.subplots(1, 1, figsize=(8,5))
        E_ax.plot(energy[:,0], energy[:,1], label='NMM', marker='o')
        E_ax.set_xlabel(r'$\beta (K^{-1})$')
        E_ax.set_ylabel(r'$E_0$')
        E_ax.annotate('N = ' + str(N) + '; g = ' + str(g), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
        E_ax.minorticks_on()
        E_fig.tight_layout()
        E_fig.savefig(os.path.join(data_path,"Energy_g" + str(g) + ".png"))
            
        g_sweep_dict_E.update({g: energy})
        energy_out.close()
        plt.close('all')
        
    # Plotting the energy versus g for varied mMax
    E_fig = plt.figure(figsize=(8,3*((len(g_sweep)+1)//2)))
    xlim = 0.
    for i, g in enumerate(g_sweep_dict_E):
        ax = E_fig.add_subplot((len(g_sweep)+1)//2, 2, i+1)
        ax.plot(g_sweep_dict_E[g][:,0], g_sweep_dict_E[g][:,1], label='NMM', marker='o')
        ax.annotate('N = ' + str(N) + '; g = ' + str(g), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
        if i == 0:
            xlim = ax.get_xlim()
        else:
            ax.set_xlim(xlim)
        ax.set_ylabel(r'$E_0$')
        ax.minorticks_on()
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
        if i+1 == 2*((len(g_sweep)+1)//2) or i+1 == 2*((len(g_sweep)+1)//2)-1:
            ax.set_xlabel(r'$\beta\ (K^{-1})$')
    E_fig.suptitle(r'$\tau\ =\ $' + str(tau))
    E_fig.tight_layout()
    E_fig.savefig("Energy_NMM_BetaSweep_tau" + str(tau) + ".png")
    
    plt.close("all")
