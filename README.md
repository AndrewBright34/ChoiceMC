# ChoiceMC
Class implentation to simulate a system of dipolar planar rotors using PIMC, or PIGs. ChoiceMC.py contains all of the documentation and explanation required to set up and run a system of rotors.

ChoiceMC_Example.py contains a sample code for a basic implementation of the ChoiceMC class, which runs the PIGS simulation for a system of 2 rotors using 9 beads, an interaction strength of 1 and a transverse potential field. It also shows how to run the exact diagonalization and NMM methods for 2 rotors.

ChoiceMC_sweep.py contains sample code to perform parametric sweeps of the number of rotors (N), the number of beads (P), the number of Monte Carlo Steps (MC_steps) and the maximum size of the free rotor eigenstate basis. It also provides some code to plot the outputs from these sweeps, and records the amount of time that is required to creat the individual rotor density matrix, the nearest neighbour density matrix and perform the MC simulation.
