# Sample code for the ChoiceMC implementation
from ChoiceMC import ChoiceMC
PIMC = ChoiceMC(m_max=15, P=9, g=1, MC_steps=1000, N=2, PIGS=True, Nskip=100, Nequilibrate=100, t=0.25)   
# Creating the probability density matrix for each rotor
PIMC.createFreeRhoMarx()
# Creating the probability density matrix for nearest neighbour interactions
PIMC.createRhoVij()
# Performing MC integration using the replica trick and swapping/unswapping
PIMC.runMCReplica(ratioTrick = False, initialState = 'random')
# Saving plots of the histograms
PIMC.plotHisto_total('left', 'middle', 'right', 'initial')
PIMC.plotHisto_N('PIMC', 'middle')
# Example of how to access class attributes
print('S2 = ', str(PIMC.S2))

# Calculating the properties by exact diagonalization, this can only be done for N=2
PIMC.runExactDiagonalization()