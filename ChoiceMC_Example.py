# Sample code for the ChoiceMC implementation
from ChoiceMC import ChoiceMC
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