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
# Setting up the variables to sweep over, these keep tau at a constant 0.0015
tau_sweep = [0.05, 0.1]
g_sweep = [0.1, 1., 2., 4.]
P_sweep = [3, 5, 7, 9, 11, 13, 15, 19, 23, 27, 31, 35]
for tau in tau_sweep:
    # Making a folder to store the current test
    time_str = "Beta_Sweep_tau"+str(tau)+"-"+str(time.gmtime().tm_year)+'-'+str(time.gmtime().tm_mon)+'-'+str(time.gmtime().tm_mday)
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
    g_sweep_dict_E_bin = {}
    # Performing the sweep over mMax
    for ig, g in enumerate(g_sweep):
        print("------------------------------------------------")
        print("Starting g = " + str(g))
        
        # Creating arrays to store the E, V and O versus g data
        energy = np.zeros((len(T_sweep),3), float)
        energy_out = open(os.path.join(data_path,"Energy_g" + str(g) + '.dat'),'w')
        energy_bin = np.zeros((len(T_sweep),3), float)
        energy_out_bin = open(os.path.join(data_path,"EnergyBin_g" + str(g) + '.dat'),'w')
    
        for i, T in enumerate(T_sweep):
            print("------------------------------------------------")
            print("Starting T = " + str(T) + "; P = " + str(P_sweep[i]))
            # Creating a ChoiceMC object for the current iteration
            PIMC = ChoiceMC(m_max=5, P=P_sweep[i], g=g, MC_steps=100000, N=2, PIGS=True, Nskip=100, Nequilibrate=100, T=T)
            print("Beta = " + str(PIMC.beta) + "; Tau = " + str(PIMC.tau))
            # Creating the probability density matrix for each rotor
            PIMC.createFreeRhoMarx()
            # Creating the probability density matrix for nearest neighbour interactions
            PIMC.createRhoVij()
            # Performing MC integration
            PIMC.runMC()
            # Saving plots of the histograms
            PIMC.plotHisto('left', "middle", 'right')
            PIMC.plotHisto('PIMC')
            
            # Storing and saving the data from the current run
            energy[i,:] = [PIMC.beta, PIMC.E_MC, PIMC.E_stdError_MC]
            energy_out.write(str(PIMC.beta) + ' ' + str(PIMC.E_MC) + ' ' + str(PIMC.E_stdError_MC) + '\n')
            energy_bin[i,:] = [PIMC.beta, PIMC.E_MC_bin, PIMC.E_stdError_MC_bin]
            energy_out_bin.write(str(PIMC.beta) + ' ' + str(PIMC.E_MC_bin) + ' ' + str(PIMC.E_stdError_MC_bin) + '\n')
                
            # Closing the remaining open plots
            plt.close('all')
    
        # Plotting
        E_fig, E_ax = plt.subplots(1, 1, figsize=(8,5))
        E_ax.plot()
        E_ax.errorbar(energy[:,0], energy[:,1], energy[:,2], label='PIGS', fmt='.-', capsize=3)
        E_ax.set_xlabel(r'$\beta (K^{-1})$')
        E_ax.set_ylabel(r'$E_0$')
        E_ax.annotate('N = 2; g = ' + str(g), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
        E_ax.minorticks_on()
        E_ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
        E_fig.tight_layout()
        E_fig.savefig(os.path.join(data_path,"Energy_g" + str(g) + ".png"))
        
        # Plotting
        E_fig, E_ax = plt.subplots(1, 1, figsize=(8,5))
        E_ax.plot()
        E_ax.errorbar(energy_bin[:,0], energy_bin[:,1], energy_bin[:,2], label='PIGS', fmt='.-', capsize=3)
        E_ax.set_xlabel(r'$\beta (K^{-1})$')
        E_ax.set_ylabel(r'$E_0$')
        E_ax.annotate('N = 2; g = ' + str(g), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
        E_ax.minorticks_on()
        E_ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
        E_fig.tight_layout()
        E_fig.savefig(os.path.join(data_path,"EnergyBin_g" + str(g) + ".png"))
            
        g_sweep_dict_E.update({g: energy})
        energy_out.close()
        g_sweep_dict_E_bin.update({g: energy_bin})
        energy_out_bin.close()
        plt.close('all')
        
    # Plotting the energy versus g for varied mMax
    E_fig = plt.figure(figsize=(8,3*((len(g_sweep)+1)//2)))
    xlim = 0.
    for i, g in enumerate(g_sweep_dict_E):
        ax = E_fig.add_subplot((len(g_sweep)+1)//2, 2, i+1)
        ax.errorbar(g_sweep_dict_E[g][:,0], g_sweep_dict_E[g][:,1], g_sweep_dict_E[g][:,2], label='PIGS', fmt='.-', capsize=3)
        ax.annotate('N = 2; g = ' + str(g), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
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
    E_fig.savefig("Energy_gBetaSweep.png")
    
    # Plotting the energy versus g for varied mMax
    E_fig = plt.figure(figsize=(8,3*((len(g_sweep)+1)//2)))
    xlim = 0.
    for i, g in enumerate(g_sweep_dict_E_bin):
        ax = E_fig.add_subplot((len(g_sweep)+1)//2, 2, i+1)
        ax.errorbar(g_sweep_dict_E_bin[g][:,0], g_sweep_dict_E_bin[g][:,1], g_sweep_dict_E_bin[g][:,2], label='PIGS', fmt='.-', capsize=3)
        ax.annotate('N = 2; g = ' + str(g), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
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
    E_fig.savefig("EnergyBin_gBetaSweep.png")
    
    plt.close("all")
    
    
