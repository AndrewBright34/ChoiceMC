# Importing the ChoiceMC class
import sys
import os
try:
    from ChoiceMC import ChoiceMC, extrapolate_E0, extrapolate_eiej, loadResult
except ModuleNotFoundError:
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))
    from ChoiceMC import ChoiceMC, extrapolate_E0, extrapolate_eiej, loadResult
import matplotlib.pyplot as plt
import numpy as np
import time
parent_dir = os.getcwd()

# Parameters
# This temperature results in beta of 4.0, which was determined to relax the system to the ground state as long as tau>0.1
T = 0.25
N_sweep = [2, 4, 16, 64, 256, 1028]
# These are chosen to sweep the tau values, to extrapolate and get a fit PIGS value
P_sweep = [39, 25, 19, 15]
N = 2

if N == 2:
    g_sweep = np.sort(np.append(np.linspace(0.01, 4, 40),np.array([0.1, 1., 2.])))
else:
    g_sweep = np.linspace(0.01, 4, 20) 

# Making a folder to store the current test
time_str = "MC_Sweep_N"+str(N)+"-"+str(time.gmtime().tm_year)+'-'+str(time.gmtime().tm_mon)+'-'+str(time.gmtime().tm_mday)
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

# Creating dictionaries to store the results from the g sweep
g_sweep_dict_E = {}
g_sweep_dict_O = {}

# Creating arrays to store the fit E and O versus g data
arrE_fit = np.zeros((len(g_sweep),2), float)
arrO_fit = np.zeros((len(g_sweep),2), float)
arrO_SmallestTau = np.zeros((len(g_sweep),2), float)
E_fit_out = open("EnergyFit_N" + str(N) + '.dat','w')
O_fit_out = open("OrienCorrFit_N" + str(N) + '.dat','w')
O_SmallestTau_out = open("OrienCorrSmallestTau_N" + str(N) + '.dat','w')

for ig, g in enumerate(g_sweep):
    print("------------------------------------------------")
    print("Starting g = " + str(g))

    # Creating arrays to store the E and O versus g data
    arrE = np.zeros((len(P_sweep),3), float)
    arrO = np.zeros((len(P_sweep),3), float)
    E_out = open(os.path.join(data_path,"Energy_g" + str(round(g,3)) + '.dat'),'w')
    O_out = open(os.path.join(data_path,"OrienCorr_g" + str(round(g,3)) + '.dat'),'w')
    
    # Performing the sweep over mMax
    for iP, P in enumerate(P_sweep):
        print("------------------------------------------------")
        print("Starting P = " + str(P))
            
        # Creating a ChoiceMC object for the current iteration
        PIMC = ChoiceMC(m_max=5, P=P, g=g, MC_steps=100000, N=N, PIGS=True, Nskip=100, Nequilibrate=100, T=T)
        print("Beta = " + str(PIMC.beta) + "; Tau = " + str(PIMC.tau))
        # Creating the probability density matrix for each rotor
        PIMC.createFreeRhoMarx()
        # Creating the probability density matrix for nearest neighbour interactions
        PIMC.createRhoVij()
        # Performing MC integration
        PIMC.runMC()
        # Saving plots of the histograms
        PIMC.plotHisto_total('left', "middle", 'right')
        PIMC.plotHisto_total('PIMC')
        
        # Storing and saving the data from the current run
        arrE[iP,:] = [PIMC.tau, PIMC.E_MC, PIMC.E_stdError_MC]
        arrO[iP,:] = [PIMC.tau, PIMC.eiej_MC, PIMC.eiej_stdError_MC]
        E_out.write(str(PIMC.tau) + ' ' + str(PIMC.E_MC) + ' ' + str(PIMC.E_stdError_MC) + '\n')
        O_out.write(str(PIMC.tau) + ' ' + str(PIMC.eiej_MC) + ' ' + str(PIMC.eiej_stdError_MC) + '\n')
            
        # Closing the remaining open plots
        plt.close('all')
    
    # Updating the dictionaries to store the raw results
    g_sweep_dict_E.update({g: arrE})
    g_sweep_dict_O.update({g: arrO})
    E_out.close()
    O_out.close()
    
    # Extrapolating to the tau = 0 limit
    Efit = extrapolate_E0(arrE, 'quadratic')
    Ofit = extrapolate_eiej(arrO, 'quadratic')
    
    # Storing the fitted results
    arrE_fit[ig, :] = [g, Efit[1]]
    arrO_fit[ig, :] = [g, Ofit[1]]
    arrO_SmallestTau[ig, :] = [g, arrO[arrO[:,0]==np.min(arrO[:,0])][0,1]]
    E_fit_out.write(str(g) + ' ' + str(Efit[1]) + '\n')
    O_fit_out.write(str(g) + ' ' + str(Ofit[1]) + '\n')
    O_SmallestTau_out.write(str(g) + ' ' + str(arrO[arrO[:,0]==np.min(arrO[:,0])][0,1]) + '\n')
    
    # Creating an array of tau for plotting a smooth line
    tauAx = np.linspace(0, np.max(arrE[:,0]), 40)
    
    # Plotting E0 vs tau
    fig, ax = plt.subplots(1, 1, figsize=(8,5))
    ax.plot(tauAx, (-1*abs(Efit[0]*tauAx**2) + Efit[1]), color='k')
    ax.plot(0, Efit[1], marker='o', color='k', label='PIGS: Fit')
    ax.errorbar(arrE[:,0], arrE[:,1], arrE[:,2], label='PIGS', fmt='o', capsize=3)
    ax.set_xlabel(r'$\tau\ (K^{-1})$')
    ax.set_ylabel(r'$E_0$'+' Per Interaction')
    ax.set_xlim((0, ax.get_xlim()[1]))
    ax.annotate('N = ' + str(N) + '; g = ' + str(round(g,3)), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
    ax.minorticks_on()
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(data_path,"Energy_g" + str(round(g,3)) + ".png"))
    
    # Plotting O vs tau
    fig, ax = plt.subplots(1, 1, figsize=(8,5))
    ax.plot(tauAx, (abs(Ofit[0]*tauAx**2) + Ofit[1]), color='k')
    ax.plot(0, Ofit[1], marker='o', color='k', label='PIGS: Fit')
    ax.errorbar(arrO[:,0], arrO[:,1], arrO[:,2], label='PIGS', fmt='o', capsize=3)
    ax.set_xlabel(r'$\tau\ (K^{-1})$')
    ax.set_ylabel("Orientational Correlation")
    ax.set_xlim((0, ax.get_xlim()[1]))
    ax.annotate('N = ' + str(N) + '; g = ' + str(round(g,3)), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
    ax.minorticks_on()
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(data_path,"OrienCorr_g" + str(round(g,3)) + ".png"))

    plt.close('all')

# Plotting E0 versus g
fig, ax = plt.subplots(1, 1, figsize=(8,5))
ax.plot(arrE_fit[:,0], arrE_fit[:,1], label='PIGS', marker='o', color='k')
ax.set_xlabel('g')
ax.set_ylabel(r'$E_0$'+' Per Interaction')
ax.annotate('N = ' + str(N), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
ax.minorticks_on()
ax.legend()
fig.tight_layout()
fig.savefig("EnergyFit_N" + str(N) + ".png")

# Plotting O versus g
fig, ax = plt.subplots(1, 1, figsize=(8,5))
ax.plot(arrO_fit[:,0], arrO_fit[:,1], label='PIGS: Fit', marker='o', color='k')
ax.set_xlabel('g')
ax.set_ylabel("Orientational Correlation")
ax.annotate('N = ' + str(N), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='left', verticalalignment='top')
ax.minorticks_on()
ax.legend()
fig.tight_layout()
fig.savefig("OrienCorrFit_N" + str(N) + ".png")

# Plotting O versus g
fig, ax = plt.subplots(1, 1, figsize=(8,5))
ax.plot(arrO_SmallestTau[:,0], arrO_SmallestTau[:,1], label='PIGS', marker='o', color='k')
ax.set_xlabel('g')
ax.set_ylabel("Orientational Correlation")
ax.annotate('N = ' + str(N), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='left', verticalalignment='top')
ax.minorticks_on()
ax.legend()
fig.tight_layout()
fig.savefig("OrienCorrSmallestTau_N" + str(N) + ".png")

E_fit_out.close()
O_fit_out.close()
plt.close("all")
    
    
