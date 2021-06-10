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

time_str = "EntanglementSweep-"+str(time.gmtime().tm_year)+'-'+str(time.gmtime().tm_mon)+'-'+str(time.gmtime().tm_mday)+'-h'+str(time.gmtime().tm_hour)
path = os.path.join(os.getcwd(), time_str)
try:
    os.mkdir(path)
except FileExistsError:
    pass
os.chdir(path)

# Setting up the variables to sweep over
g_sweep = np.linspace(0.01, 2, 10)
N=4
    
# Creating arrays to store the S2 versus g data
entanglement = np.zeros((len(g_sweep),4), float)
entanglement_out = open("S2_N"+str(N)+'.dat','w')

for ig, g in enumerate(g_sweep):
    print("------------------------------------------------")
    print("Starting g = " + str(g))
    
    # Creating a ChoiceMC object for the current iteration
    PIMC = ChoiceMC(m_max=50, P=9, g=g, MC_steps=50000, N=N, PIGS=True, Nskip=100, Nequilibrate=100)
    # Creating the probability density matrix for each rotor
    PIMC.createFreeRhoMarx()
    # Creating the probability density matrix for nearest neighbour interactions
    PIMC.createRhoVij()
    # Performing MC integration
    PIMC.runMCReplica()
    # Saving plots of the histograms
    PIMC.plotHisto('left', "middle", 'right')
    PIMC.plotHisto('PIMC')
    
    # Storing and saving the data from the current run
    entanglement[ig,:] = [g, PIMC.S2, PIMC.S2_test, PIMC.S2_err]
    entanglement_out.write(str(g) + ' ' + str(PIMC.S2) + '\n')
    
    # Closing the remaining open plots
    plt.close('all')

# Plotting
S2_fig, S2_ax = plt.subplots(1, 1, figsize=(8,5))
S2_ax.plot(entanglement[:,0], entanglement[:,1], label='original', color='#ff7f0e')
S2_ax.plot(entanglement[:,0], entanglement[:,2], label='binning', color='#1f77b4')
S2_ax.errorbar(entanglement[:,0], entanglement[:,2], entanglement[:,3], fmt='.', capsize=4, color='#1f77b4')
S2_ax.legend()
S2_ax.minorticks_on()
S2_ax.set_xlabel('g')
S2_ax.set_ylabel(r'$S_2$')
S2_ax.annotate('N = ' + str(N), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
S2_fig.tight_layout()
S2_fig.savefig("S2_N" + str(N) + ".png")

entanglement_out.close()
plt.close('all')

