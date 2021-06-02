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

time_str = "Sweep-"+str(time.gmtime().tm_year)+'-'+str(time.gmtime().tm_mon)+'-'+str(time.gmtime().tm_mday)+'-h'+str(time.gmtime().tm_hour)
path = os.path.join(os.getcwd(), time_str)
try:
    os.mkdir(path)
except FileExistsError:
    pass
os.chdir(path)

# Setting up the variables to sweep over
g_sweep = np.linspace(0.01, 8, 10)
N=4
    
# Creating arrays to store the S2 versus g data
entanglement = np.zeros((len(g_sweep),2), float)
entanglement_out = open("S2_N"+str(N)+'.dat','w')

for ig, g in enumerate(g_sweep):
    print("------------------------------------------------")
    print("Starting g = " + str(g))
    
    # Creating a ChoiceMC object for the current iteration
    PIMC = ChoiceMC(m_max=50, P=9, g=g, MC_steps=5000, N=N, PIGS=True, Nskip=100, Nequilibrate=100)
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
    entanglement[ig,:] = [g, PIMC.S2]
    entanglement_out.write(str(g) + ' ' + str(PIMC.S2) + '\n')
    
    # Closing the remaining open plots
    plt.close('all')

# Plotting
S2_fig, S2_ax = plt.subplots(1, 1, figsize=(8,5))
S2_ax.plot(entanglement[:,0], entanglement[:,1])
S2_ax.set_xlabel('Interaction Strength')
S2_ax.set_ylabel(r'Second Renyi Entropy, $S_2$')
S2_ax.set_title('Second Renyi Entropy versus Interaction Strength for ' + str(N) + " Rotor(s)")
S2_fig.tight_layout()
S2_fig.savefig("S2_N" + str(N) + ".png")

entanglement_out.close()
plt.close('all')

