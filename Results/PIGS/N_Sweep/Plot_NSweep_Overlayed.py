import os
import numpy as np
import matplotlib.pyplot as plt

files = os.listdir(os.getcwd())

E_files = [file for file in files if "Energy" in file and '.dat' in file]
O_files = [file for file in files if "OrienCorr" in file and '.dat' in file]

E_dict = {}
for file in E_files:
    E_dict.update({int((file.split('_N')[1]).split('.dat')[0]) : np.loadtxt(file)})

O_dict = {}
for file in O_files:
    O_dict.update({int((file.split('_N')[1]).split('.dat')[0]) : np.loadtxt(file)})

# Plotting E0 versus g
fig, ax = plt.subplots(1, 1, figsize=(8,5))
for N in sorted(E_dict):
    ax.plot(E_dict[N][:,0], E_dict[N][:,1], label='N = ' + str(N), marker='o')
ax.set_xlabel('g')
ax.set_ylabel(r'$E_0$'+' Per Interaction')
ax.minorticks_on()
ax.legend()
fig.tight_layout()
fig.savefig("EnergyFitOverall.png")

# Plotting O versus g
fig, ax = plt.subplots(1, 1, figsize=(8,5))
for N in sorted(O_dict):
    ax.plot(O_dict[N][:,0], O_dict[N][:,1], label='N = ' + str(N), marker='o')
ax.set_xlabel('g')
ax.set_ylabel("Orientational Correlation")
ax.minorticks_on()
ax.legend()
fig.tight_layout()
fig.savefig("OrienCorrFitOverall.png")

plt.close('all')