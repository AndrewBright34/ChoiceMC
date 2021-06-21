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
time_str = "ED_Sweep-"+str(time.gmtime().tm_year)+'-'+str(time.gmtime().tm_mon)+'-'+str(time.gmtime().tm_mday)
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
g_sweep = np.sort(np.append(np.linspace(0.01, 4, 40),np.array([0.1, 1., 2.])))
mMax_sweep = [1, 2, 5, 10, 15, 25]

# Creating dictionaries to store the results from the g and N sweep
mMax_sweep_dict_E = {}
mMax_sweep_dict_V = {}
mMax_sweep_dict_O = {}

# Performing the sweep over mMax
for imMax, mMax in enumerate(mMax_sweep):
    print("------------------------------------------------")
    print("Starting m_max = " + str(mMax))
    
    # Creating arrays to store the E, V and O versus g data
    energy = np.zeros((len(g_sweep),2), float)
    energy_out = open(os.path.join(data_path,"Energy_mMax" + str(mMax) + '.dat'),'w')
    correlations = np.zeros((len(g_sweep),2), float)
    correlations_out = open(os.path.join(data_path,"OrienCorr_mMax" + str(mMax) + '.dat'),'w')
    for ig, g in enumerate(g_sweep):
        print("------------------------------------------------")
        print("Starting g = " + str(g))
        # Creating a ChoiceMC object for the current iteration
        PIMC = ChoiceMC(m_max=mMax, P=9, g=g, MC_steps=10000, N=2, PIGS=True, Nskip=100, Nequilibrate=100)
        PIMC.exactDiagonalization()
        
        # Storing and saving the data from the current run
        energy[ig,:] = [g, PIMC.E0_ED]
        correlations[ig,:] = [g, PIMC.eiej_ED]
        energy_out.write(str(g) + ' ' + str(PIMC.E0_ED) + '\n')
        correlations_out.write(str(g) + ' ' + str(PIMC.eiej_ED) + '\n')
            
        # Closing the remaining open plots
        plt.close('all')
    
    # Plotting
    E_fig, E_ax = plt.subplots(1, 1, figsize=(8,5))
    E_ax.plot(energy[:,0], energy[:,1])
    E_ax.set_xlabel('g')
    E_ax.set_ylabel(r'$E_0$')
    E_ax.annotate('N = 2; ' + r'$M_{Max} = $' + str(mMax), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
    E_ax.minorticks_on()
    E_ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    E_fig.tight_layout()
    E_fig.savefig(os.path.join(data_path,"Energy_mMax" + str(mMax) + ".png"))
    
    O_fig, O_ax = plt.subplots(1, 1, figsize=(8,5))
    O_ax.plot(correlations[:,0], correlations[:,1])
    O_ax.set_xlabel('g')
    O_ax.set_ylabel('Orientational Correlation')
    O_ax.annotate('N = 2; ' + r'$M_{Max} = $' + str(mMax), xy=(0.5, 0.10),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
    O_ax.minorticks_on()
    O_ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    O_fig.tight_layout()
    O_fig.savefig(os.path.join(data_path,"OrienCorr_mMax" + str(mMax) + ".png"))
        
    mMax_sweep_dict_E.update({mMax: energy})
    mMax_sweep_dict_O.update({mMax: correlations})

    energy_out.close()
    correlations_out.close()
    plt.close('all')
    
# Plotting the energy versus g for varied mMax
E_fig = plt.figure(figsize=(8,3*((len(mMax_sweep)+1)//2)))
xlim = 0.
for i, mMax in enumerate(mMax_sweep_dict_E):
    ax = E_fig.add_subplot((len(mMax_sweep)+1)//2, 2, i+1)
    ax.plot(mMax_sweep_dict_E[mMax][:,0], mMax_sweep_dict_E[mMax][:,1])
    ax.annotate('N = 2; ' + r'$M_{Max} = $' + str(mMax), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
    if i == 0:
        xlim = ax.get_xlim()
    else:
        ax.set_xlim(xlim)
    ax.set_ylabel(r'$E_0$')
    ax.minorticks_on()
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    if i+1 == 2*((len(mMax_sweep)+1)//2) or i+1 == 2*((len(mMax_sweep)+1)//2)-1:
        ax.set_xlabel('g')
E_fig.tight_layout()
E_fig.savefig("Energy_GSweep.png")

# Plotting the orientational correlation versus g for varied mMax
O_fig = plt.figure(figsize=(8,3*((len(mMax_sweep)+1)//2)))
xlim  = 0.
for i, mMax in enumerate(mMax_sweep_dict_E):
    ax = O_fig.add_subplot((len(mMax_sweep)+1)//2, 2, i+1)
    ax.plot(mMax_sweep_dict_O[mMax][:,0], mMax_sweep_dict_O[mMax][:,1])
    ax.annotate('N = 2; ' + r'$M_{Max} = $' + str(mMax), xy=(0.5, 0.10),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
    if i == 0:
        xlim = ax.get_xlim()
    else:
        ax.set_xlim(xlim)
    ax.set_ylabel("Orientational Correlation")
    ax.minorticks_on()
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))   
    if i+1 == 2*((len(mMax_sweep)+1)//2) or i+1 == 2*((len(mMax_sweep)+1)//2)-1:
        ax.set_xlabel('g')    
O_fig.tight_layout()
O_fig.savefig("OrienCorr_GSweep.png")


# Reorganizing data to plot E0 versus mMax for convergence
g_vals = [0.1, 1., 2., 4.]
g_dict_E = {}
g_dict_O = {}
for g in g_vals:
    mMax_arr_E = np.zeros((len(mMax_sweep), 2), float)
    mMax_arr_O = np.zeros((len(mMax_sweep), 2), float)
    for i, mMax in enumerate(mMax_sweep):
        mMax_arr_E[i,:] = [mMax, mMax_sweep_dict_E[mMax][mMax_sweep_dict_E[1][:,0]==g][0,1]]
        mMax_arr_O[i,:] = [mMax, mMax_sweep_dict_O[mMax][mMax_sweep_dict_O[1][:,0]==g][0,1]]
    g_dict_E.update({g:mMax_arr_E})
    g_dict_O.update({g:mMax_arr_O})

# Plotting the energy versus g for varied mMax
E_fig = plt.figure(figsize=(8,3*((len(g_vals)+1)//2)))
xlim = 0.
for i, g in enumerate(g_dict_E):
    ax = E_fig.add_subplot((len(g_vals)+1)//2, 2, i+1)
    ax.plot(g_dict_E[g][:,0], g_dict_E[g][:,1])
    ax.annotate('N = 2; g = ' + str(g), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
    if i == 0:
        xlim = ax.get_xlim()
    else:
        ax.set_xlim(xlim)
    ax.set_ylabel(r'$E_0$')
    ax.minorticks_on()
    if i+1 == 2*((len(g_vals)+1)//2) or i+1 == 2*((len(g_vals)+1)//2)-1:
        ax.set_xlabel(r'$M_{Max}$')
E_fig.tight_layout()
E_fig.savefig("Energy_mMaxSweep.png")

# Plotting the orientational correlation versus g for varied mMax
O_fig = plt.figure(figsize=(8,3*((len(g_vals)+1)//2)))
xlim  = 0.
for i, g in enumerate(g_dict_O):
    ax = O_fig.add_subplot((len(g_vals)+1)//2, 2, i+1)
    ax.plot(g_dict_O[g][:,0], g_dict_O[g][:,1])
    ax.annotate('N = 2; g = ' + str(g), xy=(0.5, 0.10),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
    if i == 0:
        xlim = ax.get_xlim()
    else:
        ax.set_xlim(xlim)
    ax.set_ylabel("Orientational Correlation")
    ax.minorticks_on()  
    if i+1 == 2*((len(g_vals)+1)//2) or i+1 == 2*((len(g_vals)+1)//2)-1:
        ax.set_xlabel(r'$M_{Max}$')    
O_fig.tight_layout()
O_fig.savefig("OrienCorr_mMaxSweep.png")

plt.close("all")
    
    
