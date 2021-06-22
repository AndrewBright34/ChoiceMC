# Importing the ChoiceMC class
import sys
import os
try:
    from ChoiceMC import ChoiceMC, extrapolate_E0, extrapolate_eiej, loadResult
except ModuleNotFoundError:
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))
    from ChoiceMC import ChoiceMC, extrapolate_E0, extrapolate_eiej, loadResult
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import time

parent_dir = os.getcwd()
# Making a folder to store the current test
time_str = "NMM_Sweep-"+str(time.gmtime().tm_year)+'-'+str(time.gmtime().tm_mon)+'-'+str(time.gmtime().tm_mday)
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
mMax_sweep = [5, 10, 15, 25]
g_sweep = [0.1, 1.0, 2.0, 4.0]
P_sweep = range(3,42,2)

for mMax in mMax_sweep:
    g_sweep_dict = {}
    for ig, g in enumerate(g_sweep):
        energy = np.zeros((len(P_sweep),2), float)
        energy_out = open(os.path.join(data_path,"EnergyPSweep_mMax" + str(mMax) + "_g" + str(g) + '.dat'),'w')
        for iP, P in enumerate(P_sweep):
            PIMC = ChoiceMC(m_max=mMax, P=P, g=g, MC_steps=10000, N=2, PIGS=True, Nskip=100, Nequilibrate=100, T=0.5)
            PIMC.runNMM()
            energy[iP,:] = [P, PIMC.E0_NMM]
            energy_out.write(str(P) + ' ' + str(PIMC.E0_NMM) + '\n')
        
        g_sweep_dict.update({g: energy})
        
        # Plotting
        E_fig, E_ax = plt.subplots(1, 1, figsize=(8,5))
        E_ax.plot(energy[:,0], energy[:,1])
        E_ax.set_xlabel('P')
        E_ax.set_ylabel(r'$E_0$')
        E_ax.annotate('N = 2; g = ' + str(g) + r'$;\ M_{Max}$' + ' = ' + str(mMax), xy=(0.5, 0.15),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
        E_ax.minorticks_on()
        E_fig.tight_layout()
        E_fig.savefig(os.path.join(data_path,"EnergyPSweep_mMax" + str(mMax) + "_g" + str(g) + ".png"))
        
        energy_out.close()
        plt.close('all')
    
    # Plotting the energy versus g for varied mMax
    E_fig = plt.figure(figsize=(8,3*((len(g_sweep)+1)//2)))
    for i, g in enumerate(g_sweep_dict):
        ax = E_fig.add_subplot((len(g_sweep)+1)//2, 2, i+1)
        ax.plot(g_sweep_dict[g][:,0], g_sweep_dict[g][:,1])
        ax.annotate('N = 2; g = ' + str(g) + r'$;\ M_{Max}$' + ' = ' + str(mMax), xy=(0.5, 0.15),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
        ax.set_xlabel('P')
        ax.set_ylabel(r'$E_0$')
        ax.minorticks_on()
    E_fig.tight_layout()
    E_fig.savefig("Energy_PSweep_mMax" + str(mMax) + ".png")
    plt.close('all')

# Setting up the variables to sweep over
P_sweep = [3, 13, 21, 51]
g_sweep = [0.1, 1.0, 2.0, 4.0]
mMax_sweep = range(1,26)

for P in P_sweep:
    g_sweep_dict = {}
    for ig, g in enumerate(g_sweep):
        energy = np.zeros((len(mMax_sweep),2), float)
        energy_out = open(os.path.join(data_path,"EnergymMaxSweep_P" + str(P) + "_g" + str(g) + '.dat'),'w')
        for imMax, mMax in enumerate(mMax_sweep):
            PIMC = ChoiceMC(m_max=mMax, P=P, g=g, MC_steps=10000, N=2, PIGS=True, Nskip=100, Nequilibrate=100, T=0.5)
            PIMC.runNMM()
            energy[imMax,:] = [mMax, PIMC.E0_NMM]
            energy_out.write(str(mMax) + ' ' + str(PIMC.E0_NMM) + '\n')
        
        g_sweep_dict.update({g: energy})
        
        # Plotting
        E_fig, E_ax = plt.subplots(1, 1, figsize=(8,5))
        E_ax.plot(energy[:,0], energy[:,1])
        E_ax.set_xlabel(r'$M_{Max}$')
        E_ax.set_ylabel(r'$E_0$')
        E_ax.annotate('N = 2; g = ' + str(g) + '; P = 9', xy=(0.5, 0.15),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
        E_ax.minorticks_on()
        E_fig.tight_layout()
        E_fig.savefig(os.path.join(data_path,"EnergymMaxSweep_P" + str(P) + "_g" + str(g) + ".png"))
        
        energy_out.close()
        plt.close('all')
    
    # Plotting the energy versus g for varied mMax
    E_fig = plt.figure(figsize=(8,3*((len(g_sweep)+1)//2)))
    for i, g in enumerate(g_sweep_dict):
        ax = E_fig.add_subplot((len(g_sweep)+1)//2, 2, i+1)
        ax.plot(g_sweep_dict[g][:,0], g_sweep_dict[g][:,1])
        ax.annotate('N = 2; g = ' + str(g) + '; P = ' + str(P), xy=(0.5, 0.15),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
        ax.set_xlabel(r'$M_{Max}$')
        ax.set_ylabel(r'$E_0$')
        ax.minorticks_on()
    E_fig.tight_layout()
    E_fig.savefig("Energy_mMaxSweep_P" + str(P) + ".png")
    plt.close('all')

# Setting up the variables to sweep over, this is to help find the proper mMax
mMax_sweep = [5, 10]
g_sweep = [0.1, 1.0, 2.0, 4.0]
P_sweep = range(3,42,2)
T_sweep = [0.25, 0.5, 2./3., 1.]
for T in T_sweep:
    mMax_dict_E = {}
    mMax_dict_O = {}
    for mMax in mMax_sweep:
        g_sweep_dict_E = {}
        g_sweep_dict_O = {}
        for ig, g in enumerate(g_sweep):
            energy = np.zeros((len(P_sweep),2), float)
            energy_out = open(os.path.join(data_path,"Energy_mMax" + str(mMax) + "_beta" + str(round(1/T,3)) + "_g" + str(g) + '.dat'),'w')
            correlations = np.zeros((len(P_sweep),2), float)
            correlations_out = open(os.path.join(data_path,"OrienCorr_mMax" + str(mMax) + "_beta" + str(round(1/T,3)) + "_g" + str(g) + '.dat'),'w')
            for iP, P in enumerate(P_sweep):
                PIMC = ChoiceMC(m_max=mMax, P=P, g=g, MC_steps=10000, N=2, PIGS=True, Nskip=100, Nequilibrate=100, T=T)
                PIMC.runNMM()
                energy[iP,:] = [PIMC.tau, PIMC.E0_NMM]
                energy_out.write(str(PIMC.tau) + ' ' + str(PIMC.E0_NMM) + '\n')
                correlations[iP,:] = [PIMC.tau, PIMC.eiej_NMM]
                correlations_out.write(str(PIMC.tau) + ' ' + str(PIMC.eiej_NMM) + '\n')
            
            g_sweep_dict_E.update({g: energy})
            g_sweep_dict_O.update({g: correlations})
            energy_out.close()
            correlations_out.close()
            plt.close('all')
        
        # Plotting the energy versus g for varied mMax
        E_fig = plt.figure(figsize=(8,3*((len(g_sweep)+1)//2)))
        for i, g in enumerate(g_sweep_dict_E):
            ax = E_fig.add_subplot((len(g_sweep)+1)//2, 2, i+1)
            ax.plot(g_sweep_dict_E[g][:,0], g_sweep_dict_E[g][:,1])
            ax.annotate('N = 2; g = ' + str(g) + r'$;\ M_{Max}$' + ' = ' + str(mMax), xy=(0.5, 0.15),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
            ax.set_xlabel(r'$\tau\ (K^{-1})$')
            ax.set_ylabel(r'$E_0$')
            ax.minorticks_on()
        E_fig.suptitle(r'$\beta$ = ' + str(PIMC.beta) + r' $K^{-1}$')
        E_fig.tight_layout()
        E_fig.savefig("Energy_mMax" + str(mMax) + "_beta" + str(round(PIMC.beta,3)) + ".png")
        plt.close('all')
        
        # Plotting the energy versus g for varied mMax
        O_fig = plt.figure(figsize=(8,3*((len(g_sweep)+1)//2)))
        for i, g in enumerate(g_sweep_dict_O):
            ax = O_fig.add_subplot((len(g_sweep)+1)//2, 2, i+1)
            ax.plot(g_sweep_dict_O[g][:,0], g_sweep_dict_O[g][:,1])
            ax.annotate('N = 2; g = ' + str(g) + r'$;\ M_{Max}$' + ' = ' + str(mMax), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
            ax.set_xlabel(r'$\tau\ (K^{-1})$')
            ax.set_ylabel('Orientational Correlation')
            ax.minorticks_on()
        O_fig.suptitle(r'$\beta$ = ' + str(PIMC.beta) + r' $K^{-1}$')
        O_fig.tight_layout()
        O_fig.savefig("OrienCorr_mMax" + str(mMax) + "_beta" + str(round(PIMC.beta,3)) + ".png")
        plt.close('all')
        
        mMax_dict_E.update({mMax: g_sweep_dict_E})
        mMax_dict_O.update({mMax: g_sweep_dict_O})
    
    # Plotting the energy versus g for varied mMax
    E_fig = plt.figure(figsize=(8,3*((len(g_sweep)+1)//2)))
    for i, g in enumerate(g_sweep):
        ax = E_fig.add_subplot((len(g_sweep)+1)//2, 2, i+1)
        for mMax in mMax_dict_E:
            ax.plot(mMax_dict_E[mMax][g][:,0], mMax_dict_E[mMax][g][:,1], label=r'$NMM:\ M_{Max}$='+ str(mMax))
        ax.annotate('g = ' + str(g), xy=(0.5, 0.15),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
        ax.set_xlabel(r'$\tau\ (K^{-1})$')
        ax.set_ylabel(r'$E_0$')
        ax.legend(loc=1)
        ax.minorticks_on()
    E_fig.suptitle(r'N = 2; $\beta$ = ' + str(PIMC.beta) + r' $K^{-1}$')
    E_fig.tight_layout()
    E_fig.savefig("EnergyOverlayed_beta" + str(round(PIMC.beta,3)) + ".png")
    
    # Plotting the orientational correlation versus g for varied mMax
    O_fig = plt.figure(figsize=(8,3*((len(g_sweep)+1)//2)))
    for i, g in enumerate(g_sweep):
        ax = O_fig.add_subplot((len(g_sweep)+1)//2, 2, i+1)
        for mMax in mMax_dict_O:
            ax.plot(mMax_dict_O[mMax][g][:,0], mMax_dict_O[mMax][g][:,1], label=r'$NMM:\ M_{Max}$='+ str(mMax))
        ax.annotate('g = ' + str(g), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
        ax.set_xlabel(r'$\tau\ (K^{-1})$')
        ax.set_ylabel('Orientational Correlation')
        ax.legend(loc=4)
        ax.minorticks_on()
    O_fig.suptitle(r'N = 2; $\beta$ = ' + str(PIMC.beta) + r' $K^{-1}$')
    O_fig.tight_layout()
    O_fig.savefig("OrienCorrOverlayed_beta" + str(round(PIMC.beta,3)) + ".png")
    plt.close('all')

# Performing sweeps over tau and beta
tau_sweep = [0.10, 0.15, 0.20, 0.25]
g_sweep = [0.1, 1., 2., 4.]
P_sweep = range(3,46,2)
N = 2
for tau in tau_sweep:
    # Making a folder to store the current test
    time_str = "Beta_Sweep_NMM_tau"+str(tau)+"-"+str(time.gmtime().tm_year)+'-'+str(time.gmtime().tm_mon)+'-'+str(time.gmtime().tm_mday)
    path = os.path.join(parent_dir, time_str)
    try:
        os.mkdir(path)
    except FileExistsError:
        pass
    os.chdir(path)
    
    # Making a folder to store the individual results from each data point
    data_path = os.path.join(path, "Individual_Results_tau"+str(tau))
    try:
        os.mkdir(data_path)
    except FileExistsError:
        pass

    T_sweep = []
    for P in P_sweep:
        T_sweep.append(1./(tau*(P-1)))

    # Creating dictionaries to store the results from the g and N sweep
    g_sweep_dict_E = {}
    g_sweep_dict_O = {}
    # Performing the sweep over mMax
    for ig, g in enumerate(g_sweep):
        print("------------------------------------------------")
        print("Starting g = " + str(g))
        
        # Creating arrays to store the E versus g
        energy = np.zeros((len(T_sweep),2), float)
        energy_out = open(os.path.join(data_path,"Energy_g" + str(g) + '.dat'),'w')
        correlations = np.zeros((len(T_sweep),2), float)
        correlations_out = open(os.path.join(data_path,"OrienCorr_g" + str(g) + '.dat'),'w')
        
        for i, T in enumerate(T_sweep):
            print("------------------------------------------------")
            print("Starting T = " + str(T) + "; P = " + str(P_sweep[i]))
            # Creating a ChoiceMC object for the current iteration
            PIMC = ChoiceMC(m_max=5, P=P_sweep[i], g=g, MC_steps=1000, N=N, PIGS=True, Nskip=100, Nequilibrate=100, T=T)
            print("Beta = " + str(PIMC.beta) + "; Tau = " + str(PIMC.tau))
            # Creating the probability density matrix for each rotor
            PIMC.runNMM()
            
            # Storing and saving the data from the current run
            energy[i,:] = [PIMC.beta, PIMC.E0_NMM]
            energy_out.write(str(PIMC.beta) + ' ' + str(PIMC.E0_NMM) + '\n')
            correlations[i,:] = [PIMC.beta, PIMC.eiej_NMM]
            correlations_out.write(str(PIMC.beta) + ' ' + str(PIMC.eiej_NMM) + '\n')
                
            # Closing the remaining open plots
            plt.close('all')
    
        # Plotting
        E_fig, E_ax = plt.subplots(1, 1, figsize=(8,5))
        E_ax.plot(energy[:,0], energy[:,1], label='NMM', marker='o')
        E_ax.set_xlabel(r'$\beta (K^{-1})$')
        E_ax.set_ylabel(r'$E_0$')
        E_ax.annotate('N = ' + str(N) + '; g = ' + str(g) + r'; $\ M_{Max}$= 5', xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
        E_ax.minorticks_on()
        E_fig.tight_layout()
        E_fig.savefig(os.path.join(data_path,"Energy_g" + str(g) + ".png"))
        
        O_fig, O_ax = plt.subplots(1, 1, figsize=(8,5))
        O_ax.plot(correlations[:,0], correlations[:,1], label='NMM', marker='o')
        O_ax.set_xlabel(r'$\beta (K^{-1})$')
        O_ax.set_ylabel("Orientational Correlation")
        O_ax.annotate('N = ' + str(N) + '; g = ' + str(g) + r'; $\ M_{Max}$= 5', xy=(0.5, 0.1),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
        O_ax.minorticks_on()
        O_fig.tight_layout()
        O_fig.savefig(os.path.join(data_path,"OrienCorr_g" + str(g) + ".png"))
            
        g_sweep_dict_E.update({g: energy})
        g_sweep_dict_O.update({g: correlations})
        energy_out.close()
        correlations_out.close()
        plt.close('all')
        
    # Plotting the energy versus g for varied mMax
    E_fig = plt.figure(figsize=(8,3*((len(g_sweep)+1)//2)))
    xlim = 0.
    for i, g in enumerate(g_sweep_dict_E):
        ax = E_fig.add_subplot((len(g_sweep)+1)//2, 2, i+1)
        ax.plot(g_sweep_dict_E[g][:,0], g_sweep_dict_E[g][:,1], label='NMM', marker='o')
        ax.annotate('N = ' + str(N) + '; g = ' + str(g) + r'; $\ M_{Max}$= 5', xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
        if i == 0:
            xlim = ax.get_xlim()
        else:
            ax.set_xlim(xlim)
        ax.set_ylabel(r'$E_0$')
        ax.minorticks_on()
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
        if i+1 == 2*((len(g_sweep)+1)//2) or i+1 == 2*((len(g_sweep)+1)//2)-1:
            ax.set_xlabel(r'$\beta\ (K^{-1})$')
    E_fig.suptitle(r'$\tau\ =\ $' + str(tau))
    E_fig.tight_layout()
    E_fig.savefig("Energy_NMM_BetaSweep_tau" + str(tau) + ".png")
    
    # Plotting the orientational correlations versus g for varied mMax
    O_fig = plt.figure(figsize=(8,3*((len(g_sweep)+1)//2)))
    xlim = 0.
    for i, g in enumerate(g_sweep_dict_E):
        ax = O_fig.add_subplot((len(g_sweep)+1)//2, 2, i+1)
        ax.plot(g_sweep_dict_O[g][:,0], g_sweep_dict_O[g][:,1], label='NMM', marker='o')
        ax.annotate('N = ' + str(N) + '; g = ' + str(g) + r'; $\ M_{Max}$= 5', xy=(0.5, 0.1),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
        if i == 0:
            xlim = ax.get_xlim()
        else:
            ax.set_xlim(xlim)
        ax.set_ylabel('Orientational Correlation')
        ax.minorticks_on()
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
        if i+1 == 2*((len(g_sweep)+1)//2) or i+1 == 2*((len(g_sweep)+1)//2)-1:
            ax.set_xlabel(r'$\beta\ (K^{-1})$')
    O_fig.suptitle(r'$\tau\ =\ $' + str(tau))
    O_fig.tight_layout()
    O_fig.savefig("OrienCorr_NMM_BetaSweep_tau" + str(tau) + ".png")
    
    plt.close("all")

# P sweep is made for beta=4 and mMax=5
P_sweep = [39, 25, 19, 15]
# Sweeping over g values
g_sweep = np.sort(np.append(np.linspace(0.01, 4, 40),np.array([0.1, 1., 2.])))

# Loading Exact Diagonalization Results
arrE_ED = loadResult(os.path.join('ED', 'Energy_mMax5.dat'))
arrO_ED = loadResult(os.path.join('ED', 'OrienCorr_mMax5.dat'))

# Creating dictionaries to store the results from the g sweep
g_sweep_dict_E = {}
g_sweep_dict_O = {}

# Creating arrays to store the fit E and O versus g data
arrE_fit = np.zeros((len(g_sweep),2), float)
arrO_fit = np.zeros((len(g_sweep),2), float)
arrO_SmallestTau = np.zeros((len(g_sweep),2), float)
E_fit_out = open("EnergyFit_N2.dat",'w')
O_fit_out = open("OrienCorrFit_N2.dat",'w')
O_SmallestTau_out = open("OrienCorrSmallestTau_N2.dat",'w')

for ig, g in enumerate(g_sweep):
    print("------------------------------------------------")
    print("Starting g = " + str(g))

    # Creating arrays to store the E and O versus g data
    arrE = np.zeros((len(P_sweep),2), float)
    arrO = np.zeros((len(P_sweep),2), float)
    E_out = open(os.path.join(data_path,"Energy_g" + str(round(g,3)) + '.dat'),'w')
    O_out = open(os.path.join(data_path,"OrienCorr_g" + str(round(g,3)) + '.dat'),'w')
    
    # Performing the sweep over mMax
    for iP, P in enumerate(P_sweep):
        print("------------------------------------------------")
        print("Starting P = " + str(P))
            
        # Creating a ChoiceMC object for the current iteration
        PIMC = ChoiceMC(m_max=5, P=P, g=g, MC_steps=1000, N=2, PIGS=True, Nskip=100, Nequilibrate=100, T=0.25)
        print("Beta = " + str(PIMC.beta) + "; Tau = " + str(PIMC.tau))

        # Running the NMM method
        PIMC.runNMM()
        
        # Storing and saving the data from the current run
        arrE[iP,:] = [PIMC.tau, PIMC.E0_NMM]
        arrO[iP,:] = [PIMC.tau, PIMC.eiej_NMM]
        E_out.write(str(PIMC.tau) + ' ' + str(PIMC.E0_NMM) + '\n')
        O_out.write(str(PIMC.tau) + ' ' + str(PIMC.eiej_NMM) + '\n')
            
        # Closing the remaining open plots
        plt.close('all')
    
    # Updating the dictionaries to store the raw results
    g_sweep_dict_E.update({g: arrE})
    g_sweep_dict_O.update({g: arrO})
    E_out.close()
    O_out.close()
    
    # Extrapolating to the tau = 0 limit
    #Efit = np.polyfit(arrE[:,0], arrE[:,1], 2)
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
    ax.plot(tauAx, -1*abs(Efit[0]*tauAx**2) + Efit[1], color='k')
    ax.plot(0, Efit[1], marker='o', color='k', label='NMM: Fit')
    ax.scatter(arrE[:,0], arrE[:,1], label='NMM')
    ax.plot(tauAx, np.ones(np.shape(tauAx))*arrE_ED[arrE_ED[:,0]==g,1], label = "ED", linestyle='--', color='#d62728')
    ax.set_xlabel(r'$\tau\ (K^{-1})$')
    ax.set_ylabel(r'$E_0$'+' Per Interaction')
    ax.set_xlim((0, ax.get_xlim()[1]))
    ax.annotate('N = 2; g = ' + str(round(g,3)), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
    ax.minorticks_on()
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(data_path,"Energy_g" + str(round(g,3)) + ".png"))
    
    # Plotting O vs tau
    fig, ax = plt.subplots(1, 1, figsize=(8,5))
    ax.plot(tauAx, abs(Ofit[0]*tauAx**2) + Ofit[1], color='k')
    ax.plot(0, Ofit[1], marker='o', color='k', label='NMM: Fit')
    ax.plot(tauAx, np.ones(np.shape(tauAx))*arrO_ED[arrO_ED[:,0]==g,1], label = "ED", linestyle='--', color='#d62728')
    ax.scatter(arrO[:,0], arrO[:,1], label='NMM')
    ax.set_xlabel(r'$\tau\ (K^{-1})$')
    ax.set_ylabel("Orientational Correlation")
    ax.set_xlim((0, ax.get_xlim()[1]))
    ax.annotate('N = 2; g = ' + str(round(g,3)), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
    ax.minorticks_on()
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(data_path,"OrienCorr_g" + str(round(g,3)) + ".png"))

    plt.close('all')

# Plotting E0 versus g
fig, ax = plt.subplots(1, 1, figsize=(8,5))
ax.plot(arrE_fit[:,0], arrE_fit[:,1], label='NMM', marker='o', color='k')
ax.plot(arrE_ED[:,0], arrE_ED[:,1], label='ED', marker='o', color='#d62728')
ax.set_xlabel('g')
ax.set_ylabel(r'$E_0$'+' Per Interaction')
ax.annotate('N = 2', xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
ax.minorticks_on()
ax.legend()
fig.tight_layout()
fig.savefig("Energy_N2.png")

# Plotting O versus g
fig, ax = plt.subplots(1, 1, figsize=(8,5))
ax.plot(arrO_fit[:,0], arrO_fit[:,1], label='NMM', marker='o', color='k')
ax.plot(arrO_ED[:,0], arrO_ED[:,1], label='ED', marker='o', color='#d62728')
ax.set_xlabel('g')
ax.set_ylabel("Orientational Correlation")
ax.annotate('N = 2', xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='left', verticalalignment='top')
ax.minorticks_on()
ax.legend()
fig.tight_layout()
fig.savefig("OrienCorr_N2.png")

# Plotting O versus g
fig, ax = plt.subplots(1, 1, figsize=(8,5))
ax.plot(arrO_SmallestTau[:,0], arrO_SmallestTau[:,1], label='NMM', marker='o', color='k')
ax.plot(arrO_ED[:,0], arrO_ED[:,1], label='ED', marker='o', color='#d62728')
ax.set_xlabel('g')
ax.set_ylabel("Orientational Correlation")
ax.annotate('N = 2', xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='left', verticalalignment='top')
ax.minorticks_on()
ax.legend()
fig.tight_layout()
fig.savefig("OrienCorrSmallestTau_N2.png")

E_fit_out.close()
O_fit_out.close()
O_SmallestTau_out.close()
plt.close("all")


