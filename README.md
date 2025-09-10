# AutoMMM


## 📊 Workflow Overview

The full workflow follows the scheme illustrated in the figure below:

<p align="center">
  <img width="800" alt="Workflow" src="https://github.com/user-attachments/assets/5aebb77d-33af-425a-b689-86082407c0d7" />
</p>

**Figure 2.** Workflow of our automated computational platform for constructing and analyzing MOF/polymer composites.  
It integrates **DFT calculations** for MOF slab model generation, **FF-EMD** for composite construction and gas translational/rotational dynamics, and **FF-GCMC** for gas adsorption studies.  
Interfacial exploration includes:
- **Structural characterization:** atomic density profiles, radial distribution functions, Delaunay tessellation  
- **Textural analysis:** pore size distribution, void fraction, free pore mapping  
- **Graph theory (NetworkX):** quantifies pore connectivity and identifies possible molecular transport pathways (*new implementation highlighted by the star symbol*).  


## 📂 Repository Structure
├── s1build          # Script to build the composite system

├── ALFFIVE/         # Example MOF structure folder

├── CALF20/          # Example MOF structure folder

├── Zrfcufum/        # Example MOF structure folder

├── MIL53/           # Example MOF structure folder

├── all_pols/        # Polymer libraries and related structures

├── s2gcmc           # Script to run CO₂ adsorption with CADSS (GCMC)

├── s3md             # Script to run Molecular Dynamics simulations (FF-EMD)

├── inputs_msd/      # Input files for MSD and dynamics analysis

├── s4postanalysis   # Script for post-processing (structural, textural, and graph theory analysis)

└── script/          # Utility scripts


## ⚙️ Workflow Overview

The full workflow follows the scheme illustrated in the figure:

1. **Build the composite (s1build)**  
   - Constructs MOF/polymer composites using force-field parameters (from DFT + FF-EMD).  
   - Prepares systems for adsorption and dynamics simulations.

2. **Gas adsorption (s2gcmc)**  
   - Runs **Grand Canonical Monte Carlo (GCMC)** simulations with **CADSS software**.  
   - Evaluates CO₂ adsorption isotherms in the composite structures.

3. **Molecular dynamics (s3md)**  
   - Performs **Equilibrium Molecular Dynamics (EMD)** simulations of adsorbed CO₂.  
   - Captures time-resolved structural and dynamical properties.

4. **Post-analysis (s4postanalysis)**  
   - **Structural analysis:** atomic density profiles, RDF, planarity, Delaunay tessellation.  
   - **Textural analysis:** pore size distribution, void fraction, free pore mapping.  
   - **Graph theory analysis:** connectivity, assortativity, betweenness centrality, eccentricity (via `networkx`).  
   - **Gas dynamics:** diffusivity (Green-Kubo), angular reorientation.


## 🔬 Key Features

- Automated pipeline for **CO₂ adsorption and dynamics in composite materials**.  
- Integration of **DFT and DDEC force fields** with **CADSS (GCMC)** and **LAMMPS (MD)** simulations.  
- Advanced **post-analysis toolkit** combining statistical mechanics, pore textural characterization, and graph theory.  
- Modular design – each stage (`s1`–`s4`) can be run independently or in sequence.  

## 🚀 Getting Started

1. Clone this repository:
   ```bash
   git clone https://github.com/AlejandroDiazMarquez/AutoMMM.git
   cd AutoMMM
   
2. Run the workflow step by step:
./s1build         # Build composite
./s2gcmc          # Run GCMC adsorption
./s3md            # Run MD simulations
./s4postanalysis  # Perform analysis & graph theory

3. Analysis results will be stored in their respective subdirectories.

4. 📊 Output Examples

Structural analysis: atomic density profiles, RDF, planarity.
Textural analysis: pore size distribution (PSD), void fraction, pore maps.
Graph theory: connectivity metrics (betweenness, assortativity, eccentricity).
Gas dynamics: Green-Kubo diffusivity, angular reorientation profiles.

## 🛠 Requirements

- [**CADSS**] – for Grand Canonical Monte Carlo (GCMC) adsorption simulations  
- [**LAMMPS**](https://www.lammps.org/) – for Molecular Dynamics simulations  
- **Python 3.x** with the following packages:
  - [NumPy](https://numpy.org/)  
  - [SciPy](https://scipy.org/)  
  - [Matplotlib](https://matplotlib.org/)  
  - [NetworkX](https://networkx.org/)  
  - [MDTraj](https://www.mdtraj.org/)  

---

## 📖 Citation

If you use this repository, please cite the underlying software and our work:

- Plimpton, S. (1995). *Fast Parallel Algorithms for Short-Range Molecular Dynamics.* **J. Comput. Phys. 117**(1), 1–19. https://doi.org/10.1006/jcph.1995.1039  
- McGibbon, R. T., Beauchamp, K. A., Harrigan, M. P., Klein, C., Swails, J. M., Hernández, C. X., ... & Pande, V. S. (2015). *MDTraj: A Modern Open Library for the Analysis of Molecular Dynamics Trajectories.* **Biophys. J. 109**(8), 1528–1532. https://doi.org/10.1016/j.bpj.2015.08.015  
- Hagberg, A., Swart, P., & Chult, D. S. (2008). *Exploring Network Structure, Dynamics, and Function using NetworkX.* In Proceedings of the 7th Python in Science Conference (SciPy 2008).  
- **CADSS reference:**  
  Yang, Q., & Zhong, C. (2006). *Molecular Simulation of Carbon Dioxide/Methane/Hydrogen Mixture Adsorption in Metal−Organic Frameworks.* **J. Phys. Chem. B, 110**(36), 17776–17783. https://doi.org/10.1021/jp0628753
  
If you use results or figures generated with this repository, please also cite:  
**Díaz Márquez, A. et al. (2025).**  (Submitted).  


## 👤 Author

Developed by Alejandro Díaz Márquez
Postdoctoral Researcher in Institute Charles Gerdhardh Montpellier UM5253
