# pylattice
Pylattice is a Python package for a toolkit of first-principle calculations and lattice-model calculations.

## Features
- Supplementation of a crystal structure with it's space group
- Generation of a lattice model with various shapes (e.g., chain, ladder, kagome, etc)
- A python controller for Quantum ESPRESSO and RESPACK software packages
  - Especially, RESPACK controller is under developped

## Usages
See test.py

## Installation
- Install mVMC
  - https://www.pasums.issp.u-tokyo.ac.jp/mvmc/en/
- Install Qunatum ESPRESSO
  - https://www.quantum-espresso.org

```shell
conda create -n pylattice python=3.9
conda activate pylattice
git clone https://github.com/nkitamuraQC/pylattice.git
cd pylattice
pip install -e .
```

