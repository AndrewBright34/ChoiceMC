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
entropy = np.zeros((len(g_sweep),3), float)
entropy_out = open("SecondRenyiEntropy_N"+str(N)+'.dat','w')
purity = np.zeros((len(g_sweep),3), float)
purity_out = open("Purity_N"+str(N)+'.dat','w')

if N > 2:
    entropy_RT = np.zeros((len(g_sweep),3), float)
    entropy_RT_out = open("SecondRenyiEntropy_RatioTrick_N"+str(N)+'.dat','w')
    purity_RT = np.zeros((len(g_sweep),3), float)
    purity_RT_out = open("Purity_RatioTrick_N"+str(N)+'.dat','w')
    acceptRatio_RT_dict = {}
    acceptRatio_RT_out = open("AcceptanceRatio_RatioTrick_N"+str(N)+'.dat','w')
    acceptRatioError_RT_dict = {}
    acceptRatioError_RT_out = open("AcceptanceRatioStdError_RatioTrick_N"+str(N)+'.dat','w')

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
    entropy[ig,:] = [g, PIMC.S2_MC, PIMC.S2_stdError_MC]
    entropy_out.write(str(g) + ' ' + str(PIMC.S2_MC) + ' ' + str(PIMC.S2_stdError_MC) + '\n')
    purity[ig,:] = [g, PIMC.purity_MC, PIMC.purity_stdError_MC]
    purity_out.write(str(g) + ' ' + str(PIMC.purity_MC) + ' ' + str(PIMC.purity_stdError_MC) + '\n')
    
    if N > 2:
        # Performing MC integration with the ratio trick
        PIMC.runMCReplica(ratioTrick=True)
        # Storing and saving the data from the current run
        entropy_RT[ig,:] = [g, PIMC.S2_MC, PIMC.S2_stdError_MC]
        entropy_RT_out.write(str(g) + ' ' + str(PIMC.S2_MC) + ' ' + str(PIMC.S2_stdError_MC) + '\n')
        purity_RT[ig,:] = [g, PIMC.purity_MC, PIMC.purity_stdError_MC]
        purity_RT_out.write(str(g) + ' ' + str(PIMC.purity_MC) + ' ' + str(PIMC.purity_stdError_MC) + '\n')
        acceptRatio_RT_dict.update({g: PIMC.AR_MC_arr})
        acceptRatioError_RT_dict.update({g: PIMC.AR_stdError_MC_arr})
        
        acceptRatio_RT_out.write(str(g))
        for AR in PIMC.AR_MC_arr:
            acceptRatio_RT_out.write(' ' + str(AR))
        acceptRatio_RT_out.write('\n')
        
        acceptRatioError_RT_out.write(str(g))
        for stdAR in PIMC.AR_stdError_MC_arr:
            acceptRatioError_RT_out.write( ' ' + str(stdAR))
        acceptRatioError_RT_out.write('\n')
        
    # Closing the remaining open plots
    plt.close('all')

if N==2:
    # Loading in ED results
    arrS2_ED = loadResult(os.path.join('ED', 'SecondRenyiEntropy_mMax5.dat'))
    
if N > 2:
    acceptRatio_RT = np.zeros((len(acceptRatio_RT_dict[g_sweep[0]]), len(g_sweep), 3))
    for ig, g in enumerate(acceptRatio_RT_dict):
        for iPartition, Partition in enumerate(acceptRatio_RT_dict[g]):
            acceptRatio_RT[iPartition,ig,:] = [g, acceptRatio_RT_dict[g][iPartition], acceptRatioError_RT_dict[g][iPartition]]
# Plotting
S2_fig, S2_ax = plt.subplots(1, 1, figsize=(8,5))
S2_ax.errorbar(entropy[:,0], entropy[:,1], entropy[:,2], label='PIGS', fmt='.-', capsize=3, color='k')
if N == 2:
    S2_ax.plot(arrS2_ED[:,0], arrS2_ED[:,1], label='ED', marker='o', color='#d62728')
elif N > 2:
    S2_ax.errorbar(entropy_RT[:,0], entropy_RT[:,1], entropy_RT[:,2], label='PIGS:RT', fmt='.-', capsize=3, color='#1f77b4')
S2_ax.legend()
S2_ax.minorticks_on()
S2_ax.set_xlabel('g')
S2_ax.set_ylabel(r'$S_2$')
S2_ax.annotate('N = ' + str(N), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
S2_fig.tight_layout()
S2_fig.savefig("SecondRenyiEntropy_N" + str(N) + ".png")

# Plotting
AR_fig, AR_ax = plt.subplots(1, 1, figsize=(8,5))
AR_ax.errorbar(purity[:,0], purity[:,1], purity[:,2], label='PIGS', fmt='.-', capsize=3, color='k')
if N > 2:
    AR_ax.errorbar(purity_RT[:,0], purity_RT[:,1], purity_RT[:,2], label='PIGS:RT', fmt='.-', capsize=3, color='#1f77b4')
AR_ax.minorticks_on()
AR_ax.legend()
AR_ax.set_xlabel('g')
AR_ax.set_ylabel('Purity')
AR_ax.annotate('N = ' + str(N), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
AR_fig.tight_layout()
AR_fig.savefig("Purity_N" + str(N) + ".png")

# Plotting
AR_fig, AR_ax = plt.subplots(1, 1, figsize=(8,5))
AR_ax.errorbar(purity[:,0], purity[:,1], purity[:,2], label=r'$N_{S}/N_{U}$', fmt='.-', capsize=3)
if N > 2:
    AR_ax.errorbar(purity_RT[:,0], purity_RT[:,1], purity_RT[:,2], label=r'$N_{S}/N_{U}$: RT', fmt='.-', capsize=3)
    for i in range(np.shape(acceptRatio_RT)[0]):
        AR_ax.errorbar(acceptRatio_RT[i,:,0], acceptRatio_RT[i,:,1], acceptRatio_RT[i,:,2], label=r'$N_{'+str(i+1)+'}/N_{'+str(i)+'}$: RT', fmt='.-', capsize=3)
AR_ax.minorticks_on()
AR_ax.legend()
AR_ax.set_xlabel('g')
AR_ax.set_ylabel('Acceptance Ratio')
AR_ax.annotate('N = ' + str(N), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
AR_fig.tight_layout()
AR_fig.savefig("AcceptanceRatio_N" + str(N) + ".png")

entropy_out.close()
purity_out.close()
if N > 2:
    entropy_RT_out.close()
    purity_RT_out.close()
    acceptRatio_RT_out.close()
    acceptRatioError_RT_out.close()
plt.close('all')

