# # # Sample code for the ChoiceMC implementation
# # from ChoiceMC import ChoiceMC
# PIMC = ChoiceMC(m_max=50, P=9, g=1, MC_steps=1000, N=10, PIGS=True, Nskip=100, Nequilibrate=100, V0=20, potentialField='transverse')   
# # Creating the probability density matrix for each rotor
# PIMC.createFreeRhoMarx()
# # Creating the probability density matrix for nearest neighbour interactions
# PIMC.createRhoVij()
# # Performing MC integration
# PIMC.runMC(averagePotential = True, averageEnergy = True, orientationalCorrelations = True)
# # Saving plots of the histograms
# PIMC.plotHisto('left', "middle", 'right')
# print('<E> = ', PIMC.E_MC)

from ChoiceMC import ChoiceMC
import matplotlib.pyplot as plt
import time
import numpy as np
import os
time_str = str(time.gmtime().tm_year)+'-'+str(time.gmtime().tm_mon)+'-'+str(time.gmtime().tm_hour)
path = os.path.join(os.getcwd(), time_str)
try:
    os.mkdir(path)
except FileExistsError:
    pass
os.chdir(path)

PIMC = ChoiceMC(m_max=50, P=9, g=1, MC_steps=1000, N=10, PIGS=True, Nskip=100, Nequilibrate=100, V0=20, potentialField='transverse')   
# Creating the probability density matrix for each rotor
PIMC.createFreeRhoMarx()
# Creating the probability density matrix for nearest neighbour interactions
PIMC.createRhoVij()
# Performing MC integration
PIMC.runMC(averagePotential = True, averageEnergy = True, orientationalCorrelations = True)
# Saving plots of the histograms
PIMC.plotHisto('left', "middle", 'right')
print('<E> = ', PIMC.E_MC)
# """
# Default Parameters:
#     m_max = 50
#     P = 9
#     g = 1
#     MC_steps = 5000
#     N = 2
#     PIGS = True
#     Nskip = 100
#     Nequilibrate = 100    
# """
# # Setting up the variables to sweep over
# g_sweep = np.append(np.linspace(0.01, 2, 20), np.linspace(2.1,8,10))
# N_sweep = np.linspace(2,10,10,dtype=int)
# mMax_sweep = [10, 25, 50, 100, 150, 200, 250]
# P_sweep = [3, 7, 9, 15, 31, 63]
# MCSteps_sweep = [200, 500, 1000, 5000, 10000, 15000, 20000]

# # Setting up arrays and files to store the run time for varied N, MC steps, m Max and P
# time_N = np.zeros((len(N_sweep),5), float)
# time_N_out = open('Time_N.dat','w')
# time_mMax = np.zeros((len(mMax_sweep),5), float)
# time_mMax_out = open('Time_Ngrid.dat','w')
# time_P = np.zeros((len(P_sweep),5), float)
# time_P_out = open('Time_P.dat','w')
# time_MCSteps = np.zeros((len(MCSteps_sweep),5), float)
# time_MCSteps_out = open('Time_MCSteps.dat','w')

# # Creating dictionaries to store the results from the g and N sweep
# N_sweep_dict_E = {}
# N_sweep_dict_V = {}
# N_sweep_dict_O = {}

# for iN, N in enumerate(N_sweep):
#     print("------------------------------------------------")
#     print("Starting N = " + str(N))
    
#     # Creating arrays to store the E, V and O versus g data
#     energy = np.zeros((len(g_sweep),2), float)
#     energy_out = open("Energy_N"+str(N)+'.dat','w')
#     potential = np.zeros((len(g_sweep),2), float)
#     potential_out = open("Potential_N"+str(N)+'.dat','w')
#     correlations = np.zeros((len(g_sweep),2), float)
#     correlations_out = open("O_N"+str(N)+'.dat','w')
    
#     # Variables for storing the time required for each N value
#     tFreeRho = 0
#     tRhoVij = 0
#     tRunMC = 0
    
#     for ig, g in enumerate(g_sweep):
#         print("------------------------------------------------")
#         print("Starting g = " + str(g))
#         # Creating a ChoiceMC object for the current iteration
#         PIMC = ChoiceMC(m_max=50, P=9, g=g, MC_steps=5000, N=N, PIGS=True, Nskip=100, Nequilibrate=100)
        
#         # Creating the probability density matrix for each rotor
#         t0 = time.time()
#         PIMC.createFreeRhoMarx()
#         tFreeRho = time.time()-t0
        
#         # Creating the probability density matrix for nearest neighbour interactions
#         t0 = time.time()
#         PIMC.createRhoVij()
#         tRhoVij = time.time()-t0
        
#         # Performing MC integration
#         t0 = time.time()
#         PIMC.runMC()
#         tRunMC = time.time()-t0
#         print("The previous iteration took " + str(round(tRunMC+tRhoVij+tFreeRho,3)) + " seconds")
        
#         # Saving plots of the histograms
#         PIMC.plotHisto('left', "middle", 'right')
#         PIMC.plotHisto('PIMC')
        
#         # Storing and saving the data from the current run
#         energy[ig,:] = [g, PIMC.E_MC]
#         potential[ig,:] = [g, PIMC.V_MC]
#         correlations[ig,:] = [g, PIMC.eiej_MC]
#         energy_out.write(str(g) + ' ' + str(PIMC.E_MC) + '\n')
#         potential_out.write(str(g) + ' ' + str(PIMC.V_MC) + '\n')
#         correlations_out.write(str(g) + ' ' + str(PIMC.eiej_MC) + '\n')
        
#         # Closing the remaining open plots
#         plt.close('all')
    
#     # Storing and outputting the run time
#     tTotal = tFreeRho + tRhoVij + tRunMC
#     time_N[iN,:] = [N, tFreeRho, tRhoVij, tRunMC, tTotal]
#     time_N_out.write(str(N)+' '+str(tFreeRho)+' '+str(tRhoVij)+' '+str(tRunMC)+' '+str(tTotal)+'\n')

#     # Plotting
#     E_fig, E_ax = plt.subplots(1, 1, figsize=(8,5))
#     E_ax.plot(energy[:,0], energy[:,1])
#     E_ax.set_xlabel('Interaction Strength')
#     E_ax.set_ylabel('Ground State Energy per Rotor')
#     E_ax.set_title('Ground State Energy versus Interaction Strength for ' + str(N) + " Rotor(s)")
#     E_fig.tight_layout()
#     E_fig.savefig("Energy_N" + str(N) + ".png")
    
#     V_fig, V_ax = plt.subplots(1, 1, figsize=(8,5))
#     V_ax.plot(potential[:,0], potential[:,1])
#     V_ax.set_xlabel('Interaction Strength')
#     V_ax.set_ylabel('Potential Energy per Rotor')
#     V_ax.set_title('Potential Energy versus Interaction Strength for ' + str(N) + " Rotor(s)")
#     V_fig.tight_layout()
#     V_fig.savefig("Potential_N" + str(N) + ".png")
    
#     O_fig, O_ax = plt.subplots(1, 1, figsize=(8,5))
#     O_ax.plot(correlations[:,0], correlations[:,1])
#     O_ax.set_xlabel('Interaction Strength')
#     O_ax.set_ylabel('Orientational Correlation')
#     O_ax.set_title('Orientational Correlation versus Interaction Strength for ' + str(N) + " Rotor(s)")
#     O_fig.tight_layout()
#     O_fig.savefig("O_N" + str(N) + ".png")
    
#     # Storing the results for the current N value for plotting overlapped
#     N_sweep_dict_E.update({N: energy})
#     N_sweep_dict_V.update({N: potential})
#     N_sweep_dict_O.update({N: correlations})

#     energy_out.close()
#     potential_out.close()
#     correlations_out.close()
#     plt.close('all')

# # Plotting cumulative curves
# E_fig, E_ax = plt.subplots(1, 1, figsize=(8,5))
# for N in N_sweep_dict_E:
#     E_ax.plot(N_sweep_dict_E[N][:,0], N_sweep_dict_E[N][:,1], label="N = "+str(N))
# E_ax.set_xlabel('Interaction Strength')
# E_ax.set_ylabel('Ground State Energy per Rotor')
# E_ax.set_title('Ground State Energy versus Interaction Strength')
# E_ax.legend()
# E_fig.tight_layout()
# E_fig.savefig("Energy.png")

# V_fig, V_ax = plt.subplots(1, 1, figsize=(8,5))
# for N in N_sweep_dict_V:
#     V_ax.plot(N_sweep_dict_V[N][:,0], N_sweep_dict_V[N][:,1], label="N = "+str(N))
# V_ax.set_xlabel('Interaction Strength')
# V_ax.set_ylabel('Potential Energy per Rotor')
# V_ax.set_title('Potential Energy versus Interaction Strength')
# V_ax.legend()
# V_fig.tight_layout()
# V_fig.savefig("Potential.png")

# O_fig, O_ax = plt.subplots(1, 1, figsize=(8,5))
# for N in N_sweep_dict_O:
#     O_ax.plot(N_sweep_dict_O[N][:,0], N_sweep_dict_O[N][:,1], label="N = "+str(N))
# O_ax.set_xlabel('Interaction Strength')
# O_ax.set_ylabel('Orientational Correlation per Rotor')
# O_ax.set_title('Orientational Correlation versus Interaction Strength')
# O_ax.legend()
# O_fig.tight_layout()
# O_fig.savefig("OrientationalCorrelation.png")

# plt.close("all")

# # Sweeping over the number of grid points for timing
# for iN, N in enumerate(N_sweep):
#     print("------------------------------------------------")
#     print("Starting N = " + str(N))
#     # Creating a ChoiceMC object for the current iteration
#     PIMC = ChoiceMC(m_max=50, P=9, g=1, MC_steps=5000, N=N, PIGS=True, Nskip=100, Nequilibrate=100)
#     # Creating the probability density matrix for each rotor
#     t0 = time.time()
#     PIMC.createFreeRhoMarx()
#     tFreeRho = time.time()-t0
    
#     # Creating the probability density matrix for nearest neighbour interactions
#     t0 = time.time()
#     PIMC.createRhoVij()
#     tRhoVij = time.time()-t0
    
#     # Performing MC integration
#     t0 = time.time()
#     PIMC.runMC()
#     tRunMC = time.time()-t0
#     print("The previous iteration took " + str(round(tRunMC+tRhoVij+tFreeRho,3)) + " seconds")
    
#     # Saving plots of the histograms
#     PIMC.plotHisto('left', "middle", 'right')
#     PIMC.plotHisto('PIMC')
    
#     # Closing the remaining open plots
#     plt.close('all')
    
#     # Storing and outputting the run time
#     tTotal = tFreeRho + tRhoVij + tRunMC
#     time_N[iN,:] = [N, tFreeRho, tRhoVij, tRunMC, tTotal]
#     time_N_out.write(str(N)+' '+str(tFreeRho)+' '+str(tRhoVij)+' '+str(tRunMC)+' '+str(tTotal)+'\n')
# time_N_out.close()

# # Sweeping over the maximum size of the free rotor eigenstate basis for timing
# for imMax, mMax in enumerate(mMax_sweep):
#     print("------------------------------------------------")
#     print("Starting mMax = " + str(mMax))
#     # Creating a ChoiceMC object for the current iteration
#     PIMC = ChoiceMC(m_max=mMax, P=9, g=1, MC_steps=5000, N=2, PIGS=True, Nskip=100, Nequilibrate=100)
#     # Creating the probability density matrix for each rotor
#     t0 = time.time()
#     PIMC.createFreeRhoMarx()
#     tFreeRho = time.time()-t0
    
#     # Creating the probability density matrix for nearest neighbour interactions
#     t0 = time.time()
#     PIMC.createRhoVij()
#     tRhoVij = time.time()-t0
    
#     # Performing MC integration
#     t0 = time.time()
#     PIMC.runMC()
#     tRunMC = time.time()-t0
#     print("The previous iteration took " + str(round(tRunMC+tRhoVij+tFreeRho,3)) + " seconds")
    
#     # Saving plots of the histograms
#     PIMC.plotHisto('left', "middle", 'right')
#     PIMC.plotHisto('PIMC')
    
#     # Closing the remaining open plots
#     plt.close('all')
    
#     # Storing and outputting the run time
#     tTotal = tFreeRho + tRhoVij + tRunMC
#     time_mMax[imMax,:] = [2*mMax+1, tFreeRho, tRhoVij, tRunMC, tTotal]
#     time_mMax_out.write(str(2*mMax+1)+' '+str(tFreeRho)+' '+str(tRhoVij)+' '+str(tRunMC)+' '+str(tTotal)+'\n')
# time_mMax_out.close()

# # Sweeping over the number of beads for timing
# for iP, P in enumerate(P_sweep):
#     print("------------------------------------------------")
#     print("Starting P = " + str(P))
#     # Creating a ChoiceMC object for the current iteration
#     PIMC = ChoiceMC(m_max=50, P=P, g=1, MC_steps=5000, N=2, PIGS=True, Nskip=100, Nequilibrate=100)
#     # Creating the probability density matrix for each rotor
#     t0 = time.time()
#     PIMC.createFreeRhoMarx()
#     tFreeRho = time.time()-t0
    
#     # Creating the probability density matrix for nearest neighbour interactions
#     t0 = time.time()
#     PIMC.createRhoVij()
#     tRhoVij = time.time()-t0
    
#     # Performing MC integration
#     t0 = time.time()
#     PIMC.runMC()
#     tRunMC = time.time()-t0
#     print("The previous iteration took " + str(round(tRunMC+tRhoVij+tFreeRho,3)) + " seconds")
    
#     # Saving plots of the histograms
#     PIMC.plotHisto('left', "middle", 'right')
#     PIMC.plotHisto('PIMC')
    
#     # Closing the remaining open plots
#     plt.close('all')
    
#     # Storing and outputting the run time
#     tTotal = tFreeRho + tRhoVij + tRunMC
#     time_P[iP,:] = [P, tFreeRho, tRhoVij, tRunMC, tTotal]
#     time_P_out.write(str(P)+' '+str(tFreeRho)+' '+str(tRhoVij)+' '+str(tRunMC)+' '+str(tTotal)+'\n')
# time_P_out.close()

# # Sweeping over the number Monte Carlo steps for timing
# for iMCSteps, MCSteps in enumerate(MCSteps_sweep):
#     print("------------------------------------------------")
#     print("Starting MC Steps = " + str(MCSteps))
#     # Creating a ChoiceMC object for the current iteration
#     PIMC = ChoiceMC(m_max=50, P=9, g=1, MC_steps=MCSteps, N=2, PIGS=True, Nskip=100, Nequilibrate=100)
#     # Creating the probability density matrix for each rotor
#     t0 = time.time()
#     PIMC.createFreeRhoMarx()
#     tFreeRho = time.time()-t0
    
#     # Creating the probability density matrix for nearest neighbour interactions
#     t0 = time.time()
#     PIMC.createRhoVij()
#     tRhoVij = time.time()-t0
    
#     # Performing MC integration
#     t0 = time.time()
#     PIMC.runMC()
#     tRunMC = time.time()-t0
#     print("The previous iteration took " + str(round(tRunMC+tRhoVij+tFreeRho,3)) + " seconds")
    
#     # Saving plots of the histograms
#     PIMC.plotHisto('left', "middle", 'right')
#     PIMC.plotHisto('PIMC')
    
#     # Closing the remaining open plots
#     plt.close('all')
    
#     # Storing and outputting the run time
#     tTotal = tFreeRho + tRhoVij + tRunMC
#     time_MCSteps[iMCSteps,:] = [MCSteps, tFreeRho, tRhoVij, tRunMC, tTotal]
#     time_MCSteps_out.write(str(MCSteps)+' '+str(tFreeRho)+' '+str(tRhoVij)+' '+str(tRunMC)+' '+str(tTotal)+'\n')
# time_MCSteps_out.close()



# # Plotting total time curves, this only includes the total time but the individual can be separated
# tN_fig, tN_ax = plt.subplots(1, 1, figsize=(8,5))
# tN_ax.plot(time_N[:,0], time_N[:,4])
# tN_ax.set_xlabel('Number of Rotors')
# tN_ax.set_ylabel('Time (seconds)')
# tN_ax.set_title('Total Simulation Time versus Number of Rotors')
# tN_fig.tight_layout()
# tN_fig.savefig("Time_N_total.png")

# tmMax_fig, tmMax_ax = plt.subplots(1, 1, figsize=(8,5))
# tmMax_ax.plot(time_mMax[:,0], time_mMax[:,4])
# tmMax_ax.set_xlabel('Number of Grid Points')
# tmMax_ax.set_ylabel('Time (seconds)')
# tmMax_ax.set_title('Total Simulation Time versus Number of Grid Points')
# tmMax_fig.tight_layout()
# tmMax_fig.savefig("Time_Ngrid_total.png")

# tP_fig, tP_ax = plt.subplots(1, 1, figsize=(8,5))
# tP_ax.plot(time_P[:,0], time_P[:,4])
# tP_ax.set_xlabel('Number of Beads')
# tP_ax.set_ylabel('Time (seconds)')
# tP_ax.set_title('Total Simulation Time versus Number of Beads')
# tP_fig.tight_layout()
# tP_fig.savefig("Time_P_total.png")

# tMCSteps_fig, tMCSteps_ax = plt.subplots(1, 1, figsize=(8,5))
# tMCSteps_ax.plot(time_MCSteps[:,0], time_MCSteps[:,4])
# tMCSteps_ax.set_xlabel('Number of Monte Carlo Steps')
# tMCSteps_ax.set_ylabel('Time (seconds)')
# tMCSteps_ax.set_title('Total Simulation Time versus Number of Monte Carlo Steps')
# tMCSteps_fig.tight_layout()
# tMCSteps_fig.savefig("Time_MCSteps_total.png")
# plt.close('all')



# # Plotting Free Rho time curves, this only includes the total time but the individual can be separated
# tN_fig, tN_ax = plt.subplots(1, 1, figsize=(8,5))
# tN_ax.plot(time_N[:,0], time_N[:,1])
# tN_ax.set_xlabel('Number of Rotors')
# tN_ax.set_ylabel('Time (seconds)')
# tN_ax.set_title('Individual Rotor Density Matrix Creation Time versus Number of Rotors')
# tN_fig.tight_layout()
# tN_fig.savefig("Time_N_rhoFree.png")

# tmMax_fig, tmMax_ax = plt.subplots(1, 1, figsize=(8,5))
# tmMax_ax.plot(time_mMax[:,0], time_mMax[:,1])
# tmMax_ax.set_xlabel('Number of Grid Points')
# tmMax_ax.set_ylabel('Time (seconds)')
# tmMax_ax.set_title('Individual Rotor Density Matrix Creation Time versus Number of Grid Points')
# tmMax_fig.tight_layout()
# tmMax_fig.savefig("Time_Ngrid_rhoFree.png")

# tP_fig, tP_ax = plt.subplots(1, 1, figsize=(8,5))
# tP_ax.plot(time_P[:,0], time_P[:,1])
# tP_ax.set_xlabel('Number of Beads')
# tP_ax.set_ylabel('Time (seconds)')
# tP_ax.set_title('Individual Rotor Density Matrix Creation Time versus Number of Beads')
# tP_fig.tight_layout()
# tP_fig.savefig("Time_P_rhoFree.png")

# tMCSteps_fig, tMCSteps_ax = plt.subplots(1, 1, figsize=(8,5))
# tMCSteps_ax.plot(time_MCSteps[:,0], time_MCSteps[:,1])
# tMCSteps_ax.set_xlabel('Number of Monte Carlo Steps')
# tMCSteps_ax.set_ylabel('Time (seconds)')
# tMCSteps_ax.set_title('Individual Rotor Density Matrix Creation Time versus Number of Monte Carlo Steps')
# tMCSteps_fig.tight_layout()
# tMCSteps_fig.savefig("Time_MCSteps_rhoFree.png")
# plt.close('all')



# # Plotting Free Vij time curves, this only includes the total time but the individual can be separated
# tN_fig, tN_ax = plt.subplots(1, 1, figsize=(8,5))
# tN_ax.plot(time_N[:,0], time_N[:,2])
# tN_ax.set_xlabel('Number of Rotors')
# tN_ax.set_ylabel('Time (seconds)')
# tN_ax.set_title('NN Density Matrix Creation Time versus Number of Rotors')
# tN_fig.tight_layout()
# tN_fig.savefig("Time_N_rhoVij.png")

# tmMax_fig, tmMax_ax = plt.subplots(1, 1, figsize=(8,5))
# tmMax_ax.plot(time_mMax[:,0], time_mMax[:,2])
# tmMax_ax.set_xlabel('Number of Grid Points')
# tmMax_ax.set_ylabel('Time (seconds)')
# tmMax_ax.set_title('NN Density Matrix Creation Time versus Number of Grid Points')
# tmMax_fig.tight_layout()
# tmMax_fig.savefig("Time_Ngrid_rhoVij.png")

# tP_fig, tP_ax = plt.subplots(1, 1, figsize=(8,5))
# tP_ax.plot(time_P[:,0], time_P[:,2])
# tP_ax.set_xlabel('Number of Beads')
# tP_ax.set_ylabel('Time (seconds)')
# tP_ax.set_title('NN Density Matrix Creation Time versus Number of Beads')
# tP_fig.tight_layout()
# tP_fig.savefig("Time_P_rhoVij.png")

# tMCSteps_fig, tMCSteps_ax = plt.subplots(1, 1, figsize=(8,5))
# tMCSteps_ax.plot(time_MCSteps[:,0], time_MCSteps[:,2])
# tMCSteps_ax.set_xlabel('Number of Monte Carlo Steps')
# tMCSteps_ax.set_ylabel('Time (seconds)')
# tMCSteps_ax.set_title('NN Density Matrix Creation Time versus Number of Monte Carlo Steps')
# tMCSteps_fig.tight_layout()
# tMCSteps_fig.savefig("Time_MCSteps_rhoVij.png")
# plt.close('all')



# # Plotting MC time curves, this only includes the total time but the individual can be separated
# tN_fig, tN_ax = plt.subplots(1, 1, figsize=(8,5))
# tN_ax.plot(time_N[:,0], time_N[:,3])
# tN_ax.set_xlabel('Number of Rotors')
# tN_ax.set_ylabel('Time (seconds)')
# tN_ax.set_title('Monte Carlo Integration Time versus Number of Rotors')
# tN_fig.tight_layout()
# tN_fig.savefig("Time_N_MC.png")

# tmMax_fig, tmMax_ax = plt.subplots(1, 1, figsize=(8,5))
# tmMax_ax.plot(time_mMax[:,0], time_mMax[:,3])
# tmMax_ax.set_xlabel('Number of Grid Points')
# tmMax_ax.set_ylabel('Time (seconds)')
# tmMax_ax.set_title('Monte Carlo Integration Time versus Number of Grid Points')
# tmMax_fig.tight_layout()
# tmMax_fig.savefig("Time_Ngrid_MC.png")

# tP_fig, tP_ax = plt.subplots(1, 1, figsize=(8,5))
# tP_ax.plot(time_P[:,0], time_P[:,3])
# tP_ax.set_xlabel('Number of Beads')
# tP_ax.set_ylabel('Time (seconds)')
# tP_ax.set_title('Monte Carlo Integration Time versus Number of Beads')
# tP_fig.tight_layout()
# tP_fig.savefig("Time_P_MC.png")

# tMCSteps_fig, tMCSteps_ax = plt.subplots(1, 1, figsize=(8,5))
# tMCSteps_ax.plot(time_MCSteps[:,0], time_MCSteps[:,3])
# tMCSteps_ax.set_xlabel('Number of Monte Carlo Steps')
# tMCSteps_ax.set_ylabel('Time (seconds)')
# tMCSteps_ax.set_title('Monte Carlo Integration Time versus Number of Monte Carlo Steps')
# tMCSteps_fig.tight_layout()
# tMCSteps_fig.savefig("Time_MCSteps_MC.png")
# plt.close('all')




