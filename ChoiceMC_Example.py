# Sample code for the ChoiceMC implementation
from ChoiceMC import ChoiceMC
PIMC = ChoiceMC(m_max=15, P=9, g=1, MC_steps=1000, N=2, PIGS=True, Nskip=100, Nequilibrate=100, V0=20, potentialField='transverse', t=0.25)   
# Creating the probability density matrix for each rotor
PIMC.createFreeRhoMarx()
# Creating the probability density matrix for nearest neighbour interactions
PIMC.createRhoVij()
# Performing MC integration
PIMC.runMC(averagePotential = True, averageEnergy = True, orientationalCorrelations = True, initialState = 'random')
# Saving plots of the histograms
PIMC.plotHisto_total('left', 'middle', 'right', 'initial')
PIMC.plotHisto_N('PIMC', 'middle')
# Example of how to access class attributes
print('<E> = ', PIMC.E_MC)

# Calculating the properties by exact diagonalization, this can only be done for N=2
PIMC.runExactDiagonalization()

# Calculating the properties using the NMM method, this can only be done for N=2
PIMC.runNMM()