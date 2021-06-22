# Importing the ChoiceMC class
import sys
import os
try:
    from ChoiceMC import ChoiceMC
except ModuleNotFoundError:
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))
    from ChoiceMC import ChoiceMC
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import time


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
mMax_sweep_dict_O = {}
mMax_sweep_dict_S2 = {}
mMax_sweep_dict_SvN = {}

# Performing the sweep over mMax
for imMax, mMax in enumerate(mMax_sweep):
    print("------------------------------------------------")
    print("Starting m_max = " + str(mMax))
    
    # Creating arrays to store the E, V and O versus g data
    energy = np.zeros((len(g_sweep),2), float)
    correlations = np.zeros((len(g_sweep),2), float)
    SvN = np.zeros((len(g_sweep),2), float)
    S2 = np.zeros((len(g_sweep),2), float)
    energy_out = open(os.path.join(data_path,"Energy_mMax" + str(mMax) + '.dat'),'w')
    correlations_out = open(os.path.join(data_path,"OrienCorr_mMax" + str(mMax) + '.dat'),'w')
    SvN_out = open(os.path.join(data_path,"vonNeumannEntropy_mMax" + str(mMax) + '.dat'),'w')
    S2_out = open(os.path.join(data_path,"SecondRenyiEntropy_mMax" + str(mMax) + '.dat'),'w')
    for ig, g in enumerate(g_sweep):
        print("------------------------------------------------")
        print("Starting g = " + str(g))
        # Creating a ChoiceMC object for the current iteration
        PIMC = ChoiceMC(m_max=mMax, P=9, g=g, MC_steps=10000, N=2, PIGS=True, Nskip=100, Nequilibrate=100)
        PIMC.runExactDiagonalization()
        
        # Storing and saving the data from the current run
        energy[ig,:] = [g, PIMC.E0_ED]
        correlations[ig,:] = [g, PIMC.eiej_ED]
        S2[ig,:] = [g, PIMC.S2_ED]
        SvN[ig,:] = [g, PIMC.SvN]
        energy_out.write(str(g) + ' ' + str(PIMC.E0_ED) + '\n')
        correlations_out.write(str(g) + ' ' + str(PIMC.eiej_ED) + '\n')
        S2_out.write(str(g) + ' ' + str(PIMC.S2_ED) + '\n')
        SvN_out.write(str(g) + ' ' + str(PIMC.SvN) + '\n')
    
    # Plotting
    E_fig, E_ax = plt.subplots(1, 1, figsize=(8,5))
    E_ax.plot(energy[:,0], energy[:,1])
    E_ax.set_xlabel('g')
    E_ax.set_ylabel(r'$E_0$')
    E_ax.annotate('N = 2; ' + r'$M_{Max} = $' + str(mMax), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
    E_ax.minorticks_on()
    E_fig.tight_layout()
    E_fig.savefig(os.path.join(data_path,"Energy_mMax" + str(mMax) + ".png"))
    
    O_fig, O_ax = plt.subplots(1, 1, figsize=(8,5))
    O_ax.plot(correlations[:,0], correlations[:,1])
    O_ax.set_xlabel('g')
    O_ax.set_ylabel('Orientational Correlation')
    O_ax.annotate('N = 2; ' + r'$M_{Max} = $' + str(mMax), xy=(0.5, 0.10),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
    O_ax.minorticks_on()
    O_fig.tight_layout()
    O_fig.savefig(os.path.join(data_path,"OrienCorr_mMax" + str(mMax) + ".png"))
        
    # Plotting
    S2_fig, S2_ax = plt.subplots(1, 1, figsize=(8,5))
    S2_ax.plot(S2[:,0], S2[:,1])
    S2_ax.set_xlabel('g')
    S2_ax.set_ylabel(r'$S_2$')
    S2_ax.annotate('N = 2; ' + r'$M_{Max} = $' + str(mMax), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
    S2_ax.minorticks_on()
    S2_fig.tight_layout()
    S2_fig.savefig(os.path.join(data_path,"SecondRenyiEntropy_mMax" + str(mMax) + ".png"))
    
    # Plotting
    SvN_fig, SvN_ax = plt.subplots(1, 1, figsize=(8,5))
    SvN_ax.plot(SvN[:,0], SvN[:,1])
    SvN_ax.set_xlabel('g')
    SvN_ax.set_ylabel(r'$S_{vN}$')
    SvN_ax.annotate('N = 2; ' + r'$M_{Max} = $' + str(mMax), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
    SvN_ax.minorticks_on()
    SvN_fig.tight_layout()
    SvN_fig.savefig(os.path.join(data_path,"vonNeumannEntropy_mMax" + str(mMax) + ".png"))
    
    mMax_sweep_dict_E.update({mMax: energy})
    mMax_sweep_dict_O.update({mMax: correlations})
    mMax_sweep_dict_S2.update({mMax: S2})
    mMax_sweep_dict_SvN.update({mMax: SvN})

    energy_out.close()
    correlations_out.close()
    S2_out.close()
    SvN_out.close()
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
for i, mMax in enumerate(mMax_sweep_dict_O):
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

# Plotting S2 versus g for varied mMax
S2_fig = plt.figure(figsize=(8,3*((len(mMax_sweep)+1)//2)))
xlim  = 0.
for i, mMax in enumerate(mMax_sweep_dict_S2):
    ax = S2_fig.add_subplot((len(mMax_sweep)+1)//2, 2, i+1)
    ax.plot(mMax_sweep_dict_S2[mMax][:,0], mMax_sweep_dict_S2[mMax][:,1])
    ax.annotate('N = 2; ' + r'$M_{Max} = $' + str(mMax), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
    if i == 0:
        xlim = ax.get_xlim()
    else:
        ax.set_xlim(xlim)
    ax.set_ylabel(r'$S_2$')
    ax.minorticks_on()
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))   
    if i+1 == 2*((len(mMax_sweep)+1)//2) or i+1 == 2*((len(mMax_sweep)+1)//2)-1:
        ax.set_xlabel('g')    
S2_fig.tight_layout()
S2_fig.savefig("SecondRenyiEntropy_GSweep.png")

# Plotting SvN versus g for varied mMax
SvN_fig = plt.figure(figsize=(8,3*((len(mMax_sweep)+1)//2)))
xlim  = 0.
for i, mMax in enumerate(mMax_sweep_dict_S2):
    ax = SvN_fig.add_subplot((len(mMax_sweep)+1)//2, 2, i+1)
    ax.plot(mMax_sweep_dict_SvN[mMax][:,0], mMax_sweep_dict_SvN[mMax][:,1])
    ax.annotate('N = 2; ' + r'$M_{Max} = $' + str(mMax), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
    if i == 0:
        xlim = ax.get_xlim()
    else:
        ax.set_xlim(xlim)
    ax.set_ylabel(r'$S_{vN}$')
    ax.minorticks_on()
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))   
    if i+1 == 2*((len(mMax_sweep)+1)//2) or i+1 == 2*((len(mMax_sweep)+1)//2)-1:
        ax.set_xlabel('g')    
SvN_fig.tight_layout()
SvN_fig.savefig("vonNeumannEntropy_GSweep.png")


# Reorganizing data to plot E0 versus mMax for convergence
g_vals = [0.1, 1., 2., 4.]
g_dict_E = {}
g_dict_O = {}
g_dict_S2 = {}
g_dict_SvN = {}
for g in g_vals:
    mMax_arr_E = np.zeros((len(mMax_sweep), 2), float)
    mMax_arr_O = np.zeros((len(mMax_sweep), 2), float)
    mMax_arr_S2 = np.zeros((len(mMax_sweep), 2), float)
    mMax_arr_SvN = np.zeros((len(mMax_sweep), 2), float)
    for i, mMax in enumerate(mMax_sweep):
        mMax_arr_E[i,:] = [mMax, mMax_sweep_dict_E[mMax][mMax_sweep_dict_E[1][:,0]==g][0,1]]
        mMax_arr_O[i,:] = [mMax, mMax_sweep_dict_O[mMax][mMax_sweep_dict_O[1][:,0]==g][0,1]]
        mMax_arr_S2[i,:] = [mMax, mMax_sweep_dict_S2[mMax][mMax_sweep_dict_S2[1][:,0]==g][0,1]]
        mMax_arr_SvN[i,:] = [mMax, mMax_sweep_dict_SvN[mMax][mMax_sweep_dict_SvN[1][:,0]==g][0,1]]
    g_dict_E.update({g:mMax_arr_E})
    g_dict_O.update({g:mMax_arr_O})
    g_dict_S2.update({g:mMax_arr_S2})
    g_dict_SvN.update({g:mMax_arr_SvN})

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

# Plotting S2 versus g for varied mMax
S2_fig = plt.figure(figsize=(8,3*((len(g_vals)+1)//2)))
xlim  = 0.
for i, g in enumerate(g_dict_S2):
    ax = S2_fig.add_subplot((len(g_vals)+1)//2, 2, i+1)
    ax.plot(g_dict_S2[g][:,0], g_dict_S2[g][:,1])
    ax.annotate('N = 2; g = ' + str(g), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
    if i == 0:
        xlim = ax.get_xlim()
    else:
        ax.set_xlim(xlim)
    ax.set_ylabel(r'$S_{2}$')
    ax.minorticks_on()  
    if i+1 == 2*((len(g_vals)+1)//2) or i+1 == 2*((len(g_vals)+1)//2)-1:
        ax.set_xlabel(r'$M_{Max}$')    
S2_fig.tight_layout()
S2_fig.savefig("SecondRenyiEntropy_mMaxSweep.png")

# Plotting SvN versus g for varied mMax
SvN_fig = plt.figure(figsize=(8,3*((len(g_vals)+1)//2)))
xlim  = 0.
for i, g in enumerate(g_dict_SvN):
    ax = SvN_fig.add_subplot((len(g_vals)+1)//2, 2, i+1)
    ax.plot(g_dict_SvN[g][:,0], g_dict_SvN[g][:,1])
    ax.annotate('N = 2; g = ' + str(g), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
    if i == 0:
        xlim = ax.get_xlim()
    else:
        ax.set_xlim(xlim)
    ax.set_ylabel(r'$S_{vN}$')
    ax.minorticks_on()  
    if i+1 == 2*((len(g_vals)+1)//2) or i+1 == 2*((len(g_vals)+1)//2)-1:
        ax.set_xlabel(r'$M_{Max}$')    
SvN_fig.tight_layout()
SvN_fig.savefig("vonNeumannEntropy_mMaxSweep.png")

plt.close("all")
    
    
