# Sample code for the ChoiceMC implementation
from ChoiceMC import ChoiceMC
PIMC = ChoiceMC(m_max=15, P=9, g=1, MC_steps=10000, N=2, PIGS=True, Nskip=100, Nequilibrate=100, T=0.25)   
# Creating the probability density matrix for each rotor
PIMC.createFreeRhoMarx()
# Creating the probability density matrix for nearest neighbour interactions
PIMC.createRhoVij()
# Performing MC integration using the replica trick and swapping/unswapping
PIMC.runMCReplica(ratioTrick = False, initialState = 'random')
# Performing MC integration using the replica and ratio trick and swapping/unswapping
PIMC.runMCReplica(ratioTrick = True, initialState = 'random')
# Example of how to access class attributes
print('S2 = ', str(PIMC.S2_MC))

# Calculating the properties by exact diagonalization, this can only be done for N=2
PIMC.runExactDiagonalization()