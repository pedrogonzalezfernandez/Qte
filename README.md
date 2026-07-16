<img src="Qte_Logo.png" alt="Qte Logo" width="512" />

# Qte (QuantumTimeEvolution)

Qte is a set of Max/MSP externals for building finite-dimensional quantum systems and computing their time development from a given Hamiltonian.

This repository contains:
- Source code for the externals (`source/`)
- Max help patches (`help/`)
- Compiled Max/MSP externals (`external/`)
- Contributed examples and related materials (`contrib/`)

ATTENTION: The externals will only work on Mac OS.

In order to install it, copy the "qte" folder inside the "packages" folder in your Max folder. The Max packages folder is under "Documents/Max 8" (for Max 8) or "Documents/Max 9" (for Max 9).

If you installed qte "by hand" on a Mac, chances are that you'll have to remove the quarantine manually as well. The simplest way to do it is by going in the terminal and typing:

```bash
xattr -rd com.apple.quarantine ~/Documents/Max\ 8/Packages/qte/*
