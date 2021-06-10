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
time_str = "Beta_Sweep-"+str(time.gmtime().tm_year)+'-'+str(time.gmtime().tm_mon)+'-'+str(time.gmtime().tm_mday)
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

# Setting up the variables to sweep over, these keep tau at a constant 0.0015
g_sweep = [0.1, 1., 2., 4.]
T_sweep = [222.2222222222220, 133.3333333333330, 95.2380952380952, 74.0740740740741, 60.6060606060606, 
           51.2820512820513, 44.4444444444444, 35.0877192982456, 28.9855072463768, 24.6913580246914, 
           21.5053763440860, 19.0476190476190, 17.0940170940171, 16.2601626016260, 14.8148148148148, 
           13.6054421768707, 12.5786163522013, 11.6959064327485, 10.9289617486339, 10.2564102564103]
P_sweep = [3, 5, 7, 9, 11, 13, 15, 19, 23, 27, 31, 35, 39, 41, 45, 49, 53, 57, 61, 65]

# Creating dictionaries to store the results from the g and N sweep
g_sweep_dict_E = {}

# Performing the sweep over mMax
for ig, g in enumerate(g_sweep):
    print("------------------------------------------------")
    print("Starting g = " + str(g))
    
    # Creating arrays to store the E, V and O versus g data
    energy = np.zeros((len(T_sweep),3), float)
    energy_out = open(os.path.join(data_path,"Energy_g" + str(g) + '.dat'),'w')

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
    E_ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
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
    ax.errorbar(g_sweep_dict_E[g][:,0], g_sweep_dict_E[g][:,1], g_sweep_dict_E[g][:,2], label='PIGS', fmt='.-', capsize=3)
    ax.annotate('N = 2; g = ' + str(g), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
    if i == 0:
        xlim = ax.get_xlim()
    else:
        ax.set_xlim(xlim)
    ax.set_ylabel(r'$E_0$')
    ax.minorticks_on()
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    if i+1 == 2*((len(g_sweep)+1)//2) or i+1 == 2*((len(g_sweep)+1)//2)-1:
        ax.set_xlabel(r'$\beta\ (K^{-1})$')
E_fig.tight_layout()
E_fig.savefig("Energy_gBetaSweep.png")

plt.close("all")
    
    