# Importing the ChoiceMC class
import sys
import os
try:
    from ChoiceMC import ChoiceMC, loadResult
except ModuleNotFoundError:
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))
    from ChoiceMC import ChoiceMC, loadResult
import matplotlib.pyplot as plt
import time
import numpy as np

# Setting up the variables to sweep over
g_sweep = np.linspace(0.01, 3, 20)
N=4

time_str = "EntanglementSweep_N" + str(N) + "-"+str(time.gmtime().tm_year)+'-'+str(time.gmtime().tm_mon)+'-'+str(time.gmtime().tm_mday)
path = os.path.join(os.getcwd(), time_str)
try:
    os.mkdir(path)
except FileExistsError:
    pass
os.chdir(path)

# Creating arrays to store the S2 versus g data
entanglement = np.zeros((len(g_sweep),3), float)
entanglement_out = open("S2_N"+str(N)+'.dat','w')
entanglement_Ratio = np.zeros((len(g_sweep),3), float)
entanglement_Ratio_out = open("S2_RatioTrick_N"+str(N)+'.dat','w')

for ig, g in enumerate(g_sweep):
    print("------------------------------------------------")
    print("Starting g = " + str(g))
    
    # Creating a ChoiceMC object for the current iteration
    PIMC = ChoiceMC(m_max=5, P=9, g=g, MC_steps=100000, N=N, PIGS=True, Nskip=100, Nequilibrate=100, T=0.25)
    # Creating the probability density matrix for each rotor
    PIMC.createFreeRhoMarx()
    # Creating the probability density matrix for nearest neighbour interactions
    PIMC.createRhoVij()
    # Performing MC integration
    PIMC.runMCReplica()
    # Storing and saving the data from the current run
    entanglement[ig,:] = [g, PIMC.S2, PIMC.S2_stdError]
    entanglement_out.write(str(g) + ' ' + str(PIMC.S2) + ' ' + str(PIMC.S2_stdError) + '\n')
    if N > 2:
        # Performing MC integration with the ratio trick
        PIMC.runMCReplica(ratioTrick=True)
        # Storing and saving the data from the current run
        entanglement_Ratio[ig,:] = [g, PIMC.S2, PIMC.S2_stdError]
        entanglement_Ratio_out.write(str(g) + ' ' + str(PIMC.S2) + ' ' + str(PIMC.S2_stdError) + '\n')
    
    # Closing the remaining open plots
    plt.close('all')

if N==2:
    # Loading in ED results
    arrS2_ED = loadResult(os.path.join('ED', 'SecondRenyiEntropy_mMax5.dat'))

# Plotting
S2_fig, S2_ax = plt.subplots(1, 1, figsize=(8,5))
S2_ax.errorbar(entanglement[:,0], entanglement[:,1], entanglement[:,2], label='PIGS', fmt='.-', capsize=3, color='k')
if N==2:
    S2_ax.plot(arrS2_ED[:,0], arrS2_ED[:,1], label='ED', marker='o', color='#d62728')
    S2_ax.legend()
elif N>2:
    S2_ax.errorbar(entanglement_Ratio[:,0], entanglement_Ratio[:,1], entanglement_Ratio[:,2], label='PIGS:RT', fmt='.-', capsize=3, color='#1f77b4')
    S2_ax.legend()
S2_ax.minorticks_on()
S2_ax.set_xlabel('g')
S2_ax.set_ylabel(r'$S_2$')
S2_ax.annotate('N = ' + str(N), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
S2_fig.tight_layout()
S2_fig.savefig("S2_N" + str(N) + ".png")

entanglement_out.close()
plt.close('all')

