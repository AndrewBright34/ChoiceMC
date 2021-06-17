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

# Making a folder to store the individual results from each data point
data_path = os.path.join(path, "Individual_Results")
try:
    os.mkdir(data_path)
except FileExistsError:
    pass

# Setting up the variables to sweep over
g_sweep = np.linspace(0.01, 4, 30)
N_sweep = [2, 4, 8, 16, 32, 64]
mMax_sweep = [1, 2, 5, 10, 25, 50]

# Creating dictionaries to store the results from the g and N sweep
mMax_sweep_dict_E = {}
mMax_sweep_dict_V = {}
mMax_sweep_dict_O = {}

time_arr = np.zeros((len(N_sweep)*len(mMax_sweep),3))
time_out = open(os.path.join(data_path,'Time_Ngrid.dat'),'w')
time_out.write('mMax' + ' ' + 'N' + ' ' + 'Time(s)' + '\n')
t0 = time.time()

# Performing the sweep over mMax
for imMax, mMax in enumerate(mMax_sweep):
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
            
            if ig == 0:
                t0 = time.time()
            
            # Creating a ChoiceMC object for the current iteration
            PIMC = ChoiceMC(m_max=mMax, P=9, g=g, MC_steps=10000, N=N, PIGS=True, Nskip=100, Nequilibrate=100)
            
            # Creating the probability density matrix for each rotor
            PIMC.createFreeRhoMarx()
            
            # Creating the probability density matrix for nearest neighbour interactions
            PIMC.createRhoVij()
            
            # Performing MC integration
            PIMC.runMC()
            
            if ig == 0:
                MCtime = time.time()-t0
                time_out.write(str(mMax) + ' ' + str(N) + ' ' + str(MCtime) + '\n')
                time_arr[(imMax*len(mMax_sweep)+iN),:] = [mMax, N, MCtime]
                
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
        E_ax.set_xlabel('g')
        E_ax.set_ylabel(r'$E_0$'+' Per Interaction')
        E_ax.annotate('N = ' + str(N) + r'$;M_{Max} = $' + str(round(mMax,3)), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
        E_fig.tight_layout()
        E_fig.savefig(os.path.join(data_path,"Energy_N" + str(N) + "_mMax" + str(mMax) + ".png"))
        
        V_fig, V_ax = plt.subplots(1, 1, figsize=(8,5))
        V_ax.plot(potential[:,0], potential[:,1])
        V_ax.set_xlabel('g')
        V_ax.set_ylabel('<V> Per Interaction')
        V_ax.annotate('N = ' + str(N) + r'$;M_{Max} = $' + str(round(mMax,3)), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
        V_fig.tight_layout()
        V_fig.savefig(os.path.join(data_path,"Potential_N" + str(N) + "_mMax" + str(mMax) + ".png"))
        
        O_fig, O_ax = plt.subplots(1, 1, figsize=(8,5))
        O_ax.plot(correlations[:,0], correlations[:,1])
        O_ax.set_xlabel('g')
        O_ax.set_ylabel('Orientational Correlation')
        O_ax.annotate('N = ' + str(N) + r'$;M_{Max} = $' + str(round(mMax,3)), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
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
    E_ax.set_xlabel('g')
    E_ax.set_ylabel(r'$E_0$'+' Per Interaction')
    E_ax.annotate(r'$M_{Max} = $' + str(round(mMax,3)), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
    E_ax.legend()
    E_fig.tight_layout()
    E_fig.savefig("NSweep_Energy_mMax" + str(mMax) + ".png")
    
    V_fig, V_ax = plt.subplots(1, 1, figsize=(8,5))
    for N in N_sweep_dict_V:
        V_ax.plot(N_sweep_dict_V[N][:,0], N_sweep_dict_V[N][:,1], label="N = "+str(N))
    V_ax.set_xlabel('g')
    V_ax.set_ylabel('<V> Per Interaction')
    V_ax.annotate(r'$M_{Max} = $' + str(round(mMax,3)), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
    V_ax.legend()
    V_fig.tight_layout()
    V_fig.savefig("NSweep_Potential_mMax" + str(mMax) + ".png")
    
    O_fig, O_ax = plt.subplots(1, 1, figsize=(8,5))
    for N in N_sweep_dict_O:
        O_ax.plot(N_sweep_dict_O[N][:,0], N_sweep_dict_O[N][:,1], label="N = "+str(N))
    O_ax.set_xlabel('g')
    O_ax.set_ylabel('Orientational Correlation')
    O_ax.annotate(r'$M_{Max} = $' + str(round(mMax,3)), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
    O_ax.legend()
    O_fig.tight_layout()
    O_fig.savefig("NSweep_OrientationalCorrelation_mMax" + str(mMax) + ".png")
    
    plt.close("all")
    
    mMax_sweep_dict_E.update({mMax: N_sweep_dict_E})
    mMax_sweep_dict_V.update({mMax: N_sweep_dict_V})
    mMax_sweep_dict_O.update({mMax: N_sweep_dict_O})

time_out.close()

# Plotting E, V and O versus g for a fixed number of rotors and varied mMax
for N in N_sweep:
    E_fig, E_ax = plt.subplots(1, 1, figsize=(8,5))
    for mMax in mMax_sweep_dict_E:
        E_ax.plot(mMax_sweep_dict_E[mMax][N][:,0], mMax_sweep_dict_E[mMax][N][:,1], label="mMax = "+str(mMax))
    E_ax.set_xlabel('g')
    E_ax.set_ylabel(r'$E_0$'+' Per Interaction')
    E_ax.annotate('N = ' + str(N), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
    E_ax.legend()
    E_fig.tight_layout()
    E_fig.savefig("mMaxSweep_Energy_N" + str(N) + ".png")

for N in N_sweep:
    V_fig, V_ax = plt.subplots(1, 1, figsize=(8,5))
    for mMax in mMax_sweep_dict_V:
        V_ax.plot(mMax_sweep_dict_V[mMax][N][:,0], mMax_sweep_dict_V[mMax][N][:,1], label="mMax = "+str(mMax))
    V_ax.set_xlabel('g')
    V_ax.set_ylabel('<V> Per Interaction')
    V_ax.annotate('N = ' + str(N), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
    V_ax.legend()
    V_fig.tight_layout()
    V_fig.savefig("mMaxSweep_Potential_N" + str(N) + ".png")

for N in N_sweep:
    O_fig, O_ax = plt.subplots(1, 1, figsize=(8,5))
    for mMax in mMax_sweep_dict_O:
        O_ax.plot(mMax_sweep_dict_O[mMax][N][:,0], mMax_sweep_dict_O[mMax][N][:,1], label="mMax = "+str(mMax))
    O_ax.set_xlabel('g')
    O_ax.set_ylabel('Orientational Correlation')
    O_ax.annotate('N = ' + str(N), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
    O_ax.legend()
    O_fig.tight_layout()
    O_fig.savefig("mMaxSweep_OrientationalCorrelation_N" + str(N) + ".png")

# Plotting the simulation time with varied N
numRows = int(np.ceil(len(N_sweep)/3))
t_fig = plt.figure(figsize=(10,4*numRows))
for iN, N in enumerate(N_sweep):
    t_ax = t_fig.add_subplot(numRows, 3, iN+1)
    t_ax.plot(time_arr[time_arr[:,1] == N][:,0], time_arr[time_arr[:,1] == N][:,2])
    t_ax.set_xlabel(r'$M_{Max}$')
    t_ax.set_ylabel('Time (seconds)')
    t_ax.annotate('N = ' + str(N), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
t_fig.suptitle("Simulation Time versus mMax for a Varied Number of Rotors")
t_fig.tight_layout()
t_fig.savefig("mMaxSweep_Timing.png")

# Plotting the simulation time with varied mMax
numRows = int(np.ceil(len(mMax_sweep)/3))
t_fig = plt.figure(figsize=(10,4*numRows))
for imMax, mMax in enumerate(mMax_sweep):
    t_ax = t_fig.add_subplot(numRows, 3, imMax+1)
    t_ax.plot(time_arr[time_arr[:,0] == mMax][:,1], time_arr[time_arr[:,0] == mMax][:,2])
    t_ax.set_xlabel('N')
    t_ax.set_ylabel('Time (seconds)')
    t_ax.annotate(r'$M_{Max} = $' + str(mMax), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
t_fig.suptitle("Simulation Time versus the Number of Rotors for a Varied mMax")
t_fig.tight_layout()
t_fig.savefig("NSweep_Timing.png")

# Plotting the difference in simulation time from an mMax of 5
fig, ax = plt.subplots(1, 1, figsize=(8,5))
for imMax, mMax in enumerate(mMax_sweep):
    if mMax >= 5:
        ax.plot(time_arr[time_arr[:,0] == mMax][:,1], time_arr[time_arr[:,0] == mMax][:,2]-time_arr[time_arr[:,0] == 5][:,2], label=r'$M_{Max} = $' + str(mMax))
ax.set_xlabel('N')
ax.set_ylabel('Increase in Time from ' + r'$M_{Max}$' + ' = 5 (seconds)')
ax.legend()
fig.tight_layout()
fig.savefig("NSweep_TimingDifferential.png")

plt.close("all")

