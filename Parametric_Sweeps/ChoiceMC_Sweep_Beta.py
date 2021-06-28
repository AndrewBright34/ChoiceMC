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

parent_dir = os.getcwd()

N_sweep = [2,4]

# Sweeping over a couple mMax, this will ensure that mMax=5 is acceptable
for N in N_sweep:
    # Making a folder to store the current test
    time_str = "mMax_Determination_N"+str(N)+"-"+str(time.gmtime().tm_year)+'-'+str(time.gmtime().tm_mon)+'-'+str(time.gmtime().tm_mday)
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
    # Setting up the variables to sweep over, this is to help find the proper mMax
    mMax_sweep = [5, 10]
    g_sweep = [0.1, 1.0, 2.0, 4.0]
    tau_sweep = [0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5]
    T_sweep = [0.25, 0.5, 1.]
    for T in T_sweep:
        mMax_dict_E = {}
        mMax_dict_O = {}
        P_sweep = []
        for tau in tau_sweep:
            P_sweep.append(int(np.ceil(1+1/(tau*T))//2*2+1))
        P_sweep = np.unique(np.array(P_sweep))
        for mMax in mMax_sweep:
            g_sweep_dict_E = {}
            g_sweep_dict_O = {}
            for ig, g in enumerate(g_sweep):
                energy = np.zeros((len(P_sweep),3), float)
                energy_out = open(os.path.join(data_path,"Energy_mMax" + str(mMax) + "_beta" + str(round(1/T,3)) + "_g" + str(g) + '.dat'),'w')
                correlations = np.zeros((len(P_sweep),3), float)
                correlations_out = open(os.path.join(data_path,"OrienCorr_mMax" + str(mMax) + "_beta" + str(round(1/T,3)) + "_g" + str(g) + '.dat'),'w')
                for iP, P in enumerate(P_sweep):
                    PIMC = ChoiceMC(m_max=mMax, P=P, g=g, MC_steps=100000, N=N, PIGS=True, Nskip=100, Nequilibrate=100, T=T)
                    # Creating the probability density matrix for each rotor
                    PIMC.createFreeRhoMarx()
                    # Creating the probability density matrix for nearest neighbour interactions
                    PIMC.createRhoVij()
                    # Performing MC integration
                    PIMC.runMC()
                    
                    energy[iP,:] = [PIMC.tau, PIMC.E_MC, PIMC.E_stdError_MC]
                    energy_out.write(str(PIMC.tau) + ' ' + str(PIMC.E_MC) + ' ' + str(PIMC.E_stdError_MC) + '\n')
                    correlations[iP,:] = [PIMC.tau, PIMC.eiej_MC, PIMC.eiej_stdError_MC]
                    correlations_out.write(str(PIMC.tau) + ' ' + str(PIMC.eiej_MC) + ' ' + str(PIMC.eiej_stdError_MC) + '\n')
                
                g_sweep_dict_E.update({g: energy})
                g_sweep_dict_O.update({g: correlations})
                energy_out.close()
                correlations_out.close()
                plt.close('all')
            
            # Plotting the energy versus g for varied mMax
            E_fig = plt.figure(figsize=(8,3*((len(g_sweep)+1)//2)))
            for i, g in enumerate(g_sweep_dict_E):
                ax = E_fig.add_subplot((len(g_sweep)+1)//2, 2, i+1)
                ax.errorbar(g_sweep_dict_E[g][:,0], g_sweep_dict_E[g][:,1], g_sweep_dict_E[g][:,2], fmt='.-', capsize=3)
                ax.annotate('N = ' + str(N) + '; g = ' + str(g) + r'$;\ M_{Max}$' + ' = ' + str(mMax), xy=(0.5, 0.15),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
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
                ax.errorbar(g_sweep_dict_O[g][:,0], g_sweep_dict_O[g][:,1], g_sweep_dict_O[g][:,2], fmt='.-', capsize=3)
                ax.annotate('N = ' + str(N) + '; g = ' + str(g) + r'$;\ M_{Max}$' + ' = ' + str(mMax), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
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
                ax.errorbar(mMax_dict_E[mMax][g][:,0], mMax_dict_E[mMax][g][:,1], mMax_dict_E[mMax][g][:,2], fmt='.-', capsize=3, label=r'$PIGS:\ M_{Max}$='+ str(mMax))
            ax.annotate('g = ' + str(g), xy=(0.5, 0.15),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
            ax.set_xlabel(r'$\tau\ (K^{-1})$')
            ax.set_ylabel(r'$E_0$')
            ax.legend(loc=1)
            ax.minorticks_on()
        E_fig.suptitle('N = ' + str(N) + r'; $\beta$ = ' + str(PIMC.beta) + r' $K^{-1}$')
        E_fig.tight_layout()
        E_fig.savefig("EnergyOverlayed_beta" + str(round(PIMC.beta,3)) + ".png")
        
        # Plotting the orientational correlation versus g for varied mMax
        O_fig = plt.figure(figsize=(8,3*((len(g_sweep)+1)//2)))
        for i, g in enumerate(g_sweep):
            ax = O_fig.add_subplot((len(g_sweep)+1)//2, 2, i+1)
            for mMax in mMax_dict_O:
                ax.errorbar(mMax_dict_O[mMax][g][:,0], mMax_dict_O[mMax][g][:,1], mMax_dict_O[mMax][g][:,2], fmt='.-', capsize=3, label=r'$PIGS:\ M_{Max}$='+ str(mMax))
            ax.annotate('g = ' + str(g), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
            ax.set_xlabel(r'$\tau\ (K^{-1})$')
            ax.set_ylabel('Orientational Correlation')
            ax.legend(loc=4)
            ax.minorticks_on()
        O_fig.suptitle('N = ' + str(N) + r'; $\beta$ = ' + str(PIMC.beta) + r' $K^{-1}$')
        O_fig.tight_layout()
        O_fig.savefig("OrienCorrOverlayed_beta" + str(round(PIMC.beta,3)) + ".png")
        plt.close('all')

# Sweeping over beta values, to ensure that the beta selected will relax the system to the ground state
for N in N_sweep:
    # Setting up the variables to sweep over 
    tau_sweep = [0.1, 0.15, 0.2]
    g_sweep = [0.1, 1., 2., 4.]
    P_sweep = range(3,46,2)
    for tau in tau_sweep:
        # Making a folder to store the current test
        time_str = "Beta_Sweep_N"+str(N)+"_tau"+str(tau)+"-"+str(time.gmtime().tm_year)+'-'+str(time.gmtime().tm_mon)+'-'+str(time.gmtime().tm_mday)
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
            
            # Creating arrays to store the E, V and O versus g data
            energy = np.zeros((len(T_sweep),3), float)
            energy_out = open(os.path.join(data_path,"Energy_g" + str(g) + '.dat'),'w')
            correlations = np.zeros((len(T_sweep),3), float)
            correlations_out = open(os.path.join(data_path,"OrienCorr_g" + str(g) + '.dat'),'w')
        
            for i, T in enumerate(T_sweep):
                print("------------------------------------------------")
                print("Starting T = " + str(T) + "; P = " + str(P_sweep[i]))
                # Creating a ChoiceMC object for the current iteration
                PIMC = ChoiceMC(m_max=5, P=P_sweep[i], g=g, MC_steps=100000, N=N, PIGS=True, Nskip=100, Nequilibrate=100, T=T)
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
                energy[i,:] = [PIMC.beta, PIMC.E_MC, PIMC.E_stdError_MC]
                energy_out.write(str(PIMC.beta) + ' ' + str(PIMC.E_MC) + ' ' + str(PIMC.E_stdError_MC) + '\n')
                correlations[i,:] = [PIMC.beta, PIMC.eiej_MC, PIMC.eiej_stdError_MC]
                correlations_out.write(str(PIMC.beta) + ' ' + str(PIMC.eiej_MC) + ' ' + str(PIMC.eiej_stdError_MC) +  '\n')
                    
                # Closing the remaining open plots
                plt.close('all')
        
            # Plotting
            E_fig, E_ax = plt.subplots(1, 1, figsize=(8,5))
            E_ax.errorbar(energy[:,0], energy[:,1], energy[:,2], label='PIGS', fmt='.-', capsize=3)
            E_ax.set_xlabel(r'$\beta (K^{-1})$')
            E_ax.set_ylabel(r'$E_0$'+' Per Interaction')
            E_ax.annotate('N = ' + str(N) + '; g = ' + str(g), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
            E_ax.minorticks_on()
            E_fig.tight_layout()
            E_fig.savefig(os.path.join(data_path,"Energy_g" + str(g) + ".png"))
            
            O_fig, O_ax = plt.subplots(1, 1, figsize=(8,5))
            O_ax.errorbar(correlations[:,0], correlations[:,1], correlations[:,2], label='PIGS', fmt='.-', capsize=3)
            O_ax.set_xlabel(r'$\beta (K^{-1})$')
            O_ax.set_ylabel("Orientational Correlation")
            O_ax.annotate('N = ' + str(N) + '; g = ' + str(g), xy=(0.5, 0.1),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
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
            ax.errorbar(g_sweep_dict_E[g][:,0], g_sweep_dict_E[g][:,1], g_sweep_dict_E[g][:,2], label='PIGS', fmt='.-', capsize=3)
            ax.annotate('N = ' + str(N) + '; g = ' + str(g), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
            if i == 0:
                xlim = ax.get_xlim()
            else:
                ax.set_xlim(xlim)
            ax.set_ylabel(r'$E_0$'+' Per Interaction')
            ax.minorticks_on()
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
            if i+1 == 2*((len(g_sweep)+1)//2) or i+1 == 2*((len(g_sweep)+1)//2)-1:
                ax.set_xlabel(r'$\beta\ (K^{-1})$')
        E_fig.suptitle(r'$\tau\ =\ $' + str(tau))
        E_fig.tight_layout()
        E_fig.savefig("Energy_gBetaSweep_tau" + str(tau) + ".png")
        
        # Plotting the orientational correlations versus g for varied mMax
        O_fig = plt.figure(figsize=(8,3*((len(g_sweep)+1)//2)))
        xlim = 0.
        for i, g in enumerate(g_sweep_dict_O):
            ax = O_fig.add_subplot((len(g_sweep)+1)//2, 2, i+1)
            ax.errorbar(g_sweep_dict_O[g][:,0], g_sweep_dict_O[g][:,1], g_sweep_dict_O[g][:,2], label='PIGS', fmt='.-', capsize=3)
            ax.annotate('N = ' + str(N) + '; g = ' + str(g), xy=(0.5, 0.1),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
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
        O_fig.savefig("OrienCorr_gBetaSweep_tau" + str(tau) + ".png")
        
        plt.close("all")
        
    
