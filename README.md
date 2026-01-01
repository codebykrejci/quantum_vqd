
# Quantum band structure calculation using VQD and constant measurement protocol

## Introduction
This repository contains code for "Minimum measurements quantum protocol for band structure
calculations", available on https://arxiv.org/pdf/2511.04389. The code implements the **Variational Quantum Deflation (VQD)** method and **Constant Measurement Protocol** to calculate the electronic band structure of several materials using tight-binding model. 
The calculation results are stored in an **HDF5 file**, and a visualization script allows for comparison between VQD results and classical exact diagonalization.

## Requirements and Installation

To run and use this project, you need to have **Python** installed along with the packages listed in the `requirements.txt` file.
The core dependencies are:

- **Qiskit** — quantum circuits and primitives
- **SciPy** — optimization algorithms

Below is the full list of required packages as specified in the `requirements.txt` file:

```txt
contourpy==1.3.2
cycler==0.12.1
DateTime==5.5
dill==0.4.0
fonttools==4.59.0
h5py==3.14.0
kiwisolver==1.4.8
matplotlib==3.10.3
numpy==2.3.1
packaging==25.0
pbr==6.1.1
pillow==11.3.0
pyparsing==3.2.3
python-dateutil==2.9.0.post0
pytz==2025.2
qiskit==2.1.1
qiskit-aer==0.17.2
rustworkx==0.16.0
scipy==1.16.0
six==1.17.0
stevedore==5.4.1
typing_extensions==4.14.1
zope.interface==7.2
```

## Installation

Clone the repository and install the dependencies:

```bash
git clone https://github.com/codebykrejci/quantum_vqd.git
cd quantum_vqd
pip install -r requirements.txt
```

## Usage

The project is divided into two main phases: **calculation** and **visualization**.

---

### 1. Configuration Setup (`config_CuO2.py`, `config_bg.py`, `config_Si.py`)

Before running the calculation, parameters must be set in the 'config' files. There three 'config' files, one for each matarial, namely `config_CuO2.py` for 2D square lattice with three atom basis copper–oxygen CuO₂ system, `config_bg.py` for bilayer graphene and `config_Si.py` for 3D $\text{sp}^{3}\text{s}^{*}$ silicon model. The 'config' files contain the calculation parameters classical optimizer settings such as optimizer method, maximum number of iteration number of quantum circuit execution, bootstrapping. Additionally, there are material specific parameters such as the on-site energies and hopping amplitudes, high-symmetry path in the 1st Brillouin zone along which the calculation is run. The high-symmetry path is divided into a 'classical' path and a 'quantum' path. The classical path is densely sampled and treated as continuous; it is used to compute the band structure via exact diagonalization. In contrast, the quantum path consists of only a small number of discrete k-points, reflecting the higher computational cost of the VQD calculation.


### 2. Running the Calculation (`vqd_sampler.py`, `vqd_statevector.py` )

The main computation script performs the VQD calculation for every k-point in the defined band structure path. There two types of solvers. The `vqd_sampler.py` uses 'AerSimulator' method from Qiskit to sample the quantum circuit. Additionaly, the `vqd_sampler.py` imports `cmp.py` file. This file contains all necesary functions that implements the Constant Measurement Protocl as described in the main article. Furthermore, the `vqd_statevector.py` implements the VQD calculation for TB hamiltonian, Eg. (5) from the main article, using exact state vector representation of an quantum circuit.

Run the calculation using Constant Measurement Protocol:

```bash
python vqd_sampler.py
```
or the calculation using the exact state vector simulator: 
```bash
python vqd_statevector.py
```

The results are saved to an HDF5 file named `results_YYYYMMDD-HHMMSS.h5`. To change the material one only has to change 'config' option inside the `vqd_sampler.py` or `vqd_statevector.py`. 

### 3. Visualizing the Results (`vqd_printer.py`)

The visualization script allows plotting the stored results and comparing them with the reference (classical) values.
Upon execution, a file dialog will open to select the HDF5 results file.

Run the script:

```bash
python vqd_printer.py
```

The script automatically generates four plots:

- **plot_eigenvalues** — All calculated VQD energies (colored circles) vs. exact energies (gray line).
- **plot_each_eigenvalues** — Calculated energies, with each energy state plotted in a different color.
- **plot_n_fun** — The number of function evaluations per energy state.
- **plot_calc_time** — A heatmap of the minimization time (in seconds) for each state and k-point.

Users may further modify or extend the generated plots and other actions directly within the `vqd_printer.py` script.

Examples:
The files `results_20260101-175107.h5` and `results_20260101-180525.h5` contain an example results from `vqd_sampler.py` for  copper–oxygen CuO₂ system and bilayer graphen, respectively. The results can be simply opened and visualized by running `vqd_sampler.py` script. The calculation was performed with $N_{\text{max}} = 10^3$ and $N_{\text{shots}} = 2\times 10^4$.


## License
This project is released under the MIT License.

## Citation

If you use this code in your research or project, please cite it as:

```bibtex
@misc{Krejci2025Minimum,
  title        = {Minimum measurements quantum protocol for band structure calculation},
  author       = {Michal Krejčí and Lucie Krejčí and Ijaz Ahamed Mohammad and Martin Plesch and Martin Friák},
  year         = {2025},
  eprint       = {2511.04389},
  archivePrefix= {arXiv},
  primaryClass = {quant-ph},
  url          = {https://arxiv.org/abs/2511.04389},
  doi          = {10.48550/arXiv.2511.04389}
}
```


KREJČÍ, Michal, 2025. Quantum band structure calculation using VQD and constant measurement protocol [online]. GitHub repository, https://github.com/codebykrejci/quantum_vqd.git

