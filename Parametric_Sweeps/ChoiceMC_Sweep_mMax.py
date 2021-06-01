# Importing the ChoiceMC class, this structure should be improved
import sys
try:
    from ChoiceMC import ChoiceMC
except ModuleNotFoundError:
    sys.path.append('..')
    from ChoiceMC import ChoiceMC
import matplotlib.pyplot as plt
import time
import numpy as np
import os

# Making a folder to store the current test
time_str = "mMax_Sweep-"+str(time.gmtime().tm_year)+'-'+str(time.gmtime().tm_mon)+'-'+str(time.gmtime().tm_mday)
path = os.path.join(os.getcwd(), time_str)
try:
    os.mkdir(path)
except FileExistsError:
    pass
os.chdir(path)

# Making a folder to store the individual data from each data point
data_path = os.path.join(path, "Individual_Data")
try:
    os.mkdir(data_path)
except FileExistsError:
    pass

# Setting up the variables to sweep over
g_sweep = np.append(np.linspace(0.01, 2, 20), np.linspace(2.1,8,10))
N_sweep = [2, 4, 8, 16, 32, 64]
mMax_sweep = [1, 2, 5, 10, 25, 50]

# Creating dictionaries to store the results from the g and N sweep
mMax_sweep_dict_E = {}
mMax_sweep_dict_V = {}
mMax_sweep_dict_O = {}

# Performing the sweep over mMax
for mMax in mMax_sweep:
    print("------------------------------------------------")
    print("Starting m_max = " + str(mMax))
    # Creating dictionaries to store the results from the g and N sweep
    N_sweep_dict_E = {}
    N_sweep_dict_V = {}
    N_sweep_dict_O = {}
    for iN, N in enumerate(N_sweep):
        print("------------------------------------------------")
        print("Starting N = " + str(N))
        
        # Creating arrays to store the E, V and O versus g data
        energy = np.zeros((len(g_sweep),2), float)
        energy_out = open(os.path.join(data_path,"Energy_N"+str(N) + "_Mmax" + str(mMax) + '.dat'),'w')
        potential = np.zeros((len(g_sweep),2), float)
        potential_out = open(os.path.join(data_path,"Potential_N"+str(N) + "_Mmax" + str(mMax) + '.dat'),'w')
        correlations = np.zeros((len(g_sweep),2), float)
        correlations_out = open(os.path.join(data_path,"OrientationalCorrelation_N"+str(N)+ "_Mmax" + str(mMax) + '.dat'),'w')
        
        for ig, g in enumerate(g_sweep):
            print("------------------------------------------------")
            print("Starting g = " + str(g))
            # Creating a ChoiceMC object for the current iteration
            PIMC = ChoiceMC(m_max=mMax, P=9, g=g, MC_steps=10000, N=N, PIGS=True, Nskip=100, Nequilibrate=100)
            
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
            energy[ig,:] = [g, PIMC.E_MC]
            potential[ig,:] = [g, PIMC.V_MC]
            correlations[ig,:] = [g, PIMC.eiej_MC]
            energy_out.write(str(g) + ' ' + str(PIMC.E_MC) + '\n')
            potential_out.write(str(g) + ' ' + str(PIMC.V_MC) + '\n')
            correlations_out.write(str(g) + ' ' + str(PIMC.eiej_MC) + '\n')
            
            # Closing the remaining open plots
            plt.close('all')
    
        # Plotting
        E_fig, E_ax = plt.subplots(1, 1, figsize=(8,5))
        E_ax.plot(energy[:,0], energy[:,1])
        E_ax.set_xlabel('Interaction Strength')
        E_ax.set_ylabel('Ground State Energy per Rotor')
        E_ax.set_title('Ground State Energy versus Interaction Strength for ' + str(N) + " Rotors and an mMax of " + str(mMax))
        E_fig.tight_layout()
        E_fig.savefig(os.path.join(data_path,"Energy_N" + str(N) + "_mMax" + str(mMax) + ".png"))
        
        V_fig, V_ax = plt.subplots(1, 1, figsize=(8,5))
        V_ax.plot(potential[:,0], potential[:,1])
        V_ax.set_xlabel('Interaction Strength')
        V_ax.set_ylabel('Potential Energy per Rotor')
        V_ax.set_title('Potential Energy versus Interaction Strength for ' + str(N) + " Rotors and an mMax of " + str(mMax))
        V_fig.tight_layout()
        V_fig.savefig(os.path.join(data_path,"Potential_N" + str(N) + "_mMax" + str(mMax) + ".png"))
        
        O_fig, O_ax = plt.subplots(1, 1, figsize=(8,5))
        O_ax.plot(correlations[:,0], correlations[:,1])
        O_ax.set_xlabel('Interaction Strength')
        O_ax.set_ylabel('Orientational Correlation')
        O_ax.set_title('Orientational Correlation versus Interaction Strength for ' + str(N) + " Rotors and an mMax of " + str(mMax))
        O_fig.tight_layout()
        O_fig.savefig(os.path.join(data_path,"OrientationCorrelation_N" + str(N) + "_mMax" + str(mMax) + ".png"))
        
        # Storing the results for the current N value for plotting overlapped
        N_sweep_dict_E.update({N: energy})
        N_sweep_dict_V.update({N: potential})
        N_sweep_dict_O.update({N: correlations})
    
        energy_out.close()
        potential_out.close()
        correlations_out.close()
        plt.close('all')
    
    # Plotting E, V and O versus g for a fixed mMax and varied number of rotors
    E_fig, E_ax = plt.subplots(1, 1, figsize=(8,5))
    for N in N_sweep_dict_E:
        E_ax.plot(N_sweep_dict_E[N][:,0], N_sweep_dict_E[N][:,1], label="N = "+str(N))
    E_ax.set_xlabel('Interaction Strength')
    E_ax.set_ylabel('Ground State Energy per Rotor')
    E_ax.set_title('Ground State Energy versus Interaction Strength for an mMax of ' + str(mMax))
    E_ax.legend()
    E_fig.tight_layout()
    E_fig.savefig("NSweep_Energy_mMax" + str(mMax) + ".png")
    
    V_fig, V_ax = plt.subplots(1, 1, figsize=(8,5))
    for N in N_sweep_dict_V:
        V_ax.plot(N_sweep_dict_V[N][:,0], N_sweep_dict_V[N][:,1], label="N = "+str(N))
    V_ax.set_xlabel('Interaction Strength')
    V_ax.set_ylabel('Potential Energy per Rotor')
    V_ax.set_title('Potential Energy versus Interaction Strength for an mMax of ' + str(mMax))
    V_ax.legend()
    V_fig.tight_layout()
    V_fig.savefig("NSweep_Potential_mMax" + str(mMax) + ".png")
    
    O_fig, O_ax = plt.subplots(1, 1, figsize=(8,5))
    for N in N_sweep_dict_O:
        O_ax.plot(N_sweep_dict_O[N][:,0], N_sweep_dict_O[N][:,1], label="N = "+str(N))
    O_ax.set_xlabel('Interaction Strength')
    O_ax.set_ylabel('Orientational Correlation per Rotor')
    O_ax.set_title('Orientational Correlation versus Interaction Strength for an mMax of ' + str(mMax))
    O_ax.legend()
    O_fig.tight_layout()
    O_fig.savefig("NSweep_OrientationalCorrelation_mMax" + str(mMax) + ".png")
    
    plt.close("all")
    
    mMax_sweep_dict_E.update({mMax: N_sweep_dict_E})
    mMax_sweep_dict_V.update({mMax: N_sweep_dict_V})
    mMax_sweep_dict_O.update({mMax: N_sweep_dict_O})
    
# Plotting E, V and O versus g for a fixed number of rotors and varied mMax
for N in N_sweep:
    E_fig, E_ax = plt.subplots(1, 1, figsize=(8,5))
    for mMax in mMax_sweep_dict_E:
        E_ax.plot(mMax_sweep_dict_E[mMax][N][:,0], mMax_sweep_dict_E[mMax][N][:,1], label="mMax = "+str(mMax))
    E_ax.set_xlabel('Interaction Strength')
    E_ax.set_ylabel('Ground State Energy per Rotor')
    E_ax.set_title('Ground State Energy versus Interaction Strength for ' + str(N) + " Rotors")
    E_ax.legend()
    E_fig.tight_layout()
    E_fig.savefig("mMaxSweep_Energy_N" + str(N) + ".png")

for N in N_sweep:
    V_fig, V_ax = plt.subplots(1, 1, figsize=(8,5))
    for mMax in mMax_sweep_dict_V:
        V_ax.plot(mMax_sweep_dict_V[mMax][N][:,0], mMax_sweep_dict_V[mMax][N][:,1], label="mMax = "+str(mMax))
    V_ax.set_xlabel('Interaction Strength')
    V_ax.set_ylabel('Potential Energy per Rotor')
    V_ax.set_title('Potential Energy versus Interaction Strength for ' + str(N) + " Rotors")
    V_ax.legend()
    V_fig.tight_layout()
    V_fig.savefig("mMaxSweep_Potential_N" + str(N) + ".png")

for N in N_sweep:
    O_fig, O_ax = plt.subplots(1, 1, figsize=(8,5))
    for mMax in mMax_sweep_dict_O:
        O_ax.plot(mMax_sweep_dict_O[mMax][N][:,0], mMax_sweep_dict_O[mMax][N][:,1], label="mMax = "+str(mMax))
    O_ax.set_xlabel('Interaction Strength')
    O_ax.set_ylabel('Orientational Correlation per Rotor')
    O_ax.set_title('Orientational Correlation versus Interaction Strength for ' + str(N) + " Rotors")
    O_ax.legend()
    O_fig.tight_layout()
    O_fig.savefig("mMaxSweep_OrientationalCorrelation_N" + str(N) + ".png")

plt.close("all")

