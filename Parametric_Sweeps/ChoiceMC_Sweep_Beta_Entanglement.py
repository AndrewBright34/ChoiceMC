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

N_sweep = [2]

# Sweeping over a couple mMax, this will ensure that mMax=5 is acceptable
for N in N_sweep:
    # Making a folder to store the current test
    time_str = "mMax_Determination_Entanglement_N"+str(N)+"-"+str(time.gmtime().tm_year)+'-'+str(time.gmtime().tm_mon)+'-'+str(time.gmtime().tm_mday)
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
    g_sweep = [0.1, 0.25, 0.5, 1.0]
    tau_sweep = [0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5]
    T_sweep = [0.25, 0.5, 1.]
    for T in T_sweep:
        mMax_dict_S2 = {}
        mMax_dict_AR = {}
        P_sweep = []
        for tau in tau_sweep:
            P_sweep.append(int(np.ceil(1+1/(tau*T))//2*2+1))
        P_sweep = np.unique(np.array(P_sweep))
        for mMax in mMax_sweep:
            g_sweep_dict_S2 = {}
            g_sweep_dict_AR = {}
            for ig, g in enumerate(g_sweep):
                entropy = np.zeros((len(P_sweep),3), float)
                entropy_out = open(os.path.join(data_path,"SecondRenyiEntropy_mMax" + str(mMax) + "_beta" + str(round(1/T,3)) + "_g" + str(g) + '.dat'),'w')
                acceptRatio = np.zeros((len(P_sweep),2), float)
                acceptRatio_out = open(os.path.join(data_path,"AcceptRatio_mMax" + str(mMax) + "_beta" + str(round(1/T,3)) + "_g" + str(g) + '.dat'),'w')
                for iP, P in enumerate(P_sweep):
                    PIMC = ChoiceMC(m_max=mMax, P=P, g=g, MC_steps=100000, N=N, PIGS=True, Nskip=100, Nequilibrate=100, T=T)
                    # Creating the probability density matrix for each rotor
                    PIMC.createFreeRhoMarx()
                    # Creating the probability density matrix for nearest neighbour interactions
                    PIMC.createRhoVij()
                    # Performing MC integration
                    PIMC.runMCReplica()
                    
                    entropy[iP,:] = [PIMC.tau, PIMC.S2_MC, PIMC.S2_stdError_MC]
                    entropy_out.write(str(PIMC.tau) + ' ' + str(PIMC.S2_MC) + ' ' + str(PIMC.S2_stdError_MC) + '\n')
                    acceptRatio[iP,:] = [PIMC.tau, PIMC.AR_MC]
                    acceptRatio_out.write(str(PIMC.tau) + ' ' + str(PIMC.AR_MC) + '\n')
                
                g_sweep_dict_S2.update({g: entropy})
                g_sweep_dict_AR.update({g: acceptRatio})
                entropy_out.close()
                acceptRatio_out.close()
                plt.close('all')
            
            # Plotting the second renyi entropy
            S2_fig = plt.figure(figsize=(8,3*((len(g_sweep)+1)//2)))
            for i, g in enumerate(g_sweep_dict_S2):
                ax = S2_fig.add_subplot((len(g_sweep)+1)//2, 2, i+1)
                ax.errorbar(g_sweep_dict_S2[g][:,0], g_sweep_dict_S2[g][:,1], g_sweep_dict_S2[g][:,2], fmt='.-', capsize=3, color='k')
                ax.annotate('N = ' + str(N) + '; g = ' + str(g) + r'$;\ M_{Max}$' + ' = ' + str(mMax), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
                ax.set_xlabel(r'$\tau\ (K^{-1})$')
                ax.set_ylabel(r'$S_2$')
                ax.minorticks_on()
            S2_fig.suptitle(r'$\beta$ = ' + str(PIMC.beta) + r' $K^{-1}$')
            S2_fig.tight_layout()
            S2_fig.savefig("SecondRenyiEntropy_mMax" + str(mMax) + "_beta" + str(round(PIMC.beta,3)) + ".png")
            
            # Plotting the acceptance ratio
            fig = plt.figure(figsize=(8,3*((len(g_sweep)+1)//2)))
            for i, g in enumerate(g_sweep_dict_S2):
                ax = fig.add_subplot((len(g_sweep)+1)//2, 2, i+1)
                ax.plot(g_sweep_dict_AR[g][:,0], g_sweep_dict_AR[g][:,1], marker='o', color='k')
                ax.annotate('N = ' + str(N) + '; g = ' + str(g) + r'$;\ M_{Max}$' + ' = ' + str(mMax), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
                ax.set_xlabel(r'$\tau\ (K^{-1})$')
                ax.set_ylabel('Acceptance Ratio')
                ax.minorticks_on()
            fig.suptitle(r'$\beta$ = ' + str(PIMC.beta) + r' $K^{-1}$')
            fig.tight_layout()
            fig.savefig("AcceptanceRatio_mMax" + str(mMax) + "_beta" + str(round(PIMC.beta,3)) + ".png")
            
            plt.close('all')
            mMax_dict_S2.update({mMax: g_sweep_dict_S2})
            mMax_dict_AR.update({mMax: g_sweep_dict_AR})
    
        # Plotting the second renyi entropy versus g for varied mMax
        S2_fig = plt.figure(figsize=(8,3*((len(g_sweep)+1)//2)))
        for i, g in enumerate(g_sweep):
            ax = S2_fig.add_subplot((len(g_sweep)+1)//2, 2, i+1)
            for mMax in mMax_dict_S2:
                ax.errorbar(mMax_dict_S2[mMax][g][:,0], mMax_dict_S2[mMax][g][:,1], mMax_dict_S2[mMax][g][:,2], fmt='.-', capsize=3, label=r'$PIGS:\ M_{Max}$='+ str(mMax))
            ax.annotate('g = ' + str(g), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
            ax.set_xlabel(r'$\tau\ (K^{-1})$')
            ax.set_ylabel(r'$S_2$')
            ax.legend(loc=1)
            ax.minorticks_on()
        S2_fig.suptitle('N = ' + str(N) + r'; $\beta$ = ' + str(PIMC.beta) + r' $K^{-1}$')
        S2_fig.tight_layout()
        S2_fig.savefig("SecondRenyiEntropyOverlayed_beta" + str(round(PIMC.beta,3)) + ".png")
        
        # Plotting the acceptance ratio versus g for varied mMax
        fig = plt.figure(figsize=(8,3*((len(g_sweep)+1)//2)))
        for i, g in enumerate(g_sweep):
            ax = fig.add_subplot((len(g_sweep)+1)//2, 2, i+1)
            for mMax in mMax_dict_AR:
                ax.plot(mMax_dict_AR[mMax][g][:,0], mMax_dict_AR[mMax][g][:,1], marker='o', label=r'$PIGS:\ M_{Max}$='+ str(mMax))
            ax.annotate('g = ' + str(g), xy=(0.5, 0.95),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
            ax.set_xlabel(r'$\tau\ (K^{-1})$')
            ax.set_ylabel('Acceptance Ratio')
            ax.legend(loc=4)
            ax.minorticks_on()
        fig.suptitle('N = ' + str(N) + r'; $\beta$ = ' + str(PIMC.beta) + r' $K^{-1}$')
        fig.tight_layout()
        fig.savefig("AcceptanceRatioOverlayed_beta" + str(round(PIMC.beta,3)) + ".png")
        plt.close('all')

# Sweeping over beta values, to ensure that the beta selected will relax the system to the ground state
for N in N_sweep:
    # Setting up the variables to sweep over 
    tau_sweep = [0.1, 0.15, 0.2]
    g_sweep = [0.1, 1., 2., 4.]
    P_sweep = range(3,46,2)
    for tau in tau_sweep:
        # Making a folder to store the current test
        time_str = "Beta_Sweep_Entanglement_N"+str(N)+"_tau"+str(tau)+"-"+str(time.gmtime().tm_year)+'-'+str(time.gmtime().tm_mon)+'-'+str(time.gmtime().tm_mday)
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
            T_sweep.append(1./(tau*P))
    
        # Creating dictionaries to store the results from the g and N sweep
        g_sweep_dict_S2 = {}
        g_sweep_dict_AR = {}
        # Performing the sweep over mMax
        for ig, g in enumerate(g_sweep):
            print("------------------------------------------------")
            print("Starting g = " + str(g))
            
            # Creating arrays to store the E, V and O versus g data
            entropy = np.zeros((len(T_sweep),3), float)
            entropy_out = open(os.path.join(data_path,"SecondRenyiEntropy_g" + str(g) + '.dat'),'w')
            acceptRatio = np.zeros((len(T_sweep),2), float)
            acceptRatio_out = open(os.path.join(data_path,"AcceptanceRatio_g" + str(g) + '.dat'),'w')
        
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
                PIMC.runMCReplica()
                
                # Storing and saving the data from the current run
                entropy[i,:] = [PIMC.beta, PIMC.S2_MC, PIMC.S2_stdError_MC]
                entropy_out.write(str(PIMC.beta) + ' ' + str(PIMC.S2_MC) + ' ' + str(PIMC.S2_stdError_MC) + '\n')
                acceptRatio[i,:] = [PIMC.beta, PIMC.AR_MC]
                acceptRatio_out.write(str(PIMC.beta) + ' ' + str(PIMC.AR_MC) + '\n')
                    
                # Closing the remaining open plots
                plt.close('all')
        
            # Plotting
            S2_fig, S2_ax = plt.subplots(1, 1, figsize=(8,5))
            S2_ax.errorbar(entropy[:,0], entropy[:,1], entropy[:,2], label='PIGS', fmt='.-', capsize=3)
            S2_ax.set_xlabel(r'$\beta (K^{-1})$')
            S2_ax.set_ylabel(r'$S_2$')
            S2_ax.annotate('N = ' + str(N) + '; g = ' + str(g), xy=(0.5, 0.1),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
            S2_ax.minorticks_on()
            S2_fig.tight_layout()
            S2_fig.savefig(os.path.join(data_path,"SecondRenyiEntropy_g" + str(g) + ".png"))
            
            fig, ax = plt.subplots(1, 1, figsize=(8,5))
            ax.plot(acceptRatio[:,0], acceptRatio[:,1], label='PIGS', marker='o')
            ax.set_xlabel(r'$\beta (K^{-1})$')
            ax.set_ylabel("Acceptance Ratio")
            ax.annotate('N = ' + str(N) + '; g = ' + str(g), xy=(0.5, 0.1),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
            ax.minorticks_on()
            fig.tight_layout()
            fig.savefig(os.path.join(data_path,"AcceptanceRatio_g" + str(g) + ".png"))
            
            g_sweep_dict_S2.update({g: entropy})
            g_sweep_dict_AR.update({g: acceptRatio})
            entropy_out.close()
            acceptRatio_out.close()
            plt.close('all')
            
        # Plotting the energy versus g for varied mMax
        S2_fig = plt.figure(figsize=(8,3*((len(g_sweep)+1)//2)))
        xlim = 0.
        for i, g in enumerate(g_sweep_dict_S2):
            ax = S2_fig.add_subplot((len(g_sweep)+1)//2, 2, i+1)
            ax.errorbar(g_sweep_dict_S2[g][:,0], g_sweep_dict_S2[g][:,1], g_sweep_dict_S2[g][:,2], label='PIGS', fmt='.-', capsize=3)
            ax.annotate('N = ' + str(N) + '; g = ' + str(g), xy=(0.5, 0.1),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
            if i == 0:
                xlim = ax.get_xlim()
            else:
                ax.set_xlim(xlim)
            ax.set_ylabel(r'$S_2$')
            ax.minorticks_on()
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
            if i+1 == 2*((len(g_sweep)+1)//2) or i+1 == 2*((len(g_sweep)+1)//2)-1:
                ax.set_xlabel(r'$\beta\ (K^{-1})$')
        S2_fig.suptitle(r'$\tau\ =\ $' + str(tau))
        S2_fig.tight_layout()
        S2_fig.savefig("SecondRenyiEntropy_gBetaSweep" + str(tau) + ".png")
        
        # Plotting the orientational correlations versus g for varied mMax
        fig = plt.figure(figsize=(8,3*((len(g_sweep)+1)//2)))
        xlim = 0.
        for i, g in enumerate(g_sweep_dict_AR):
            ax = fig.add_subplot((len(g_sweep)+1)//2, 2, i+1)
            ax.plot(g_sweep_dict_AR[g][:,0], g_sweep_dict_AR[g][:,1], label='PIGS', marker='o')
            ax.annotate('N = ' + str(N) + '; g = ' + str(g), xy=(0.5, 0.1),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='top')
            if i == 0:
                xlim = ax.get_xlim()
            else:
                ax.set_xlim(xlim)
            ax.set_ylabel('Acceptance Ratio')
            ax.minorticks_on()
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
            if i+1 == 2*((len(g_sweep)+1)//2) or i+1 == 2*((len(g_sweep)+1)//2)-1:
                ax.set_xlabel(r'$\beta\ (K^{-1})$')
        fig.suptitle(r'$\tau\ =\ $' + str(tau))
        fig.tight_layout()
        fig.savefig("AcceptanceRatio_gBetaSweep_tau" + str(tau) + ".png")
        
        plt.close("all")
        
    
