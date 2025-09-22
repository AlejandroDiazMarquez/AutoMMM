# Practical README.txt file with basic usage instructions
MOF/Polymer Composites: CO₂ Adsorption & Dynamics Workflow
==========================================================

This tool provides a step-by-step pipeline to:
1. Build MOF/polymer composites
2. Simulate CO₂ adsorption
3. Run Molecular Dynamics
4. Perform post-analysis (structural, textural, and graph theory)

----------------------------------------------------------
Basic Instructions
----------------------------------------------------------

1. Preparation
   - Install CADSS for GCMC adsorption simulations
   - Install LAMMPS for Molecular Dynamics
   - Install Python 3.x with the following packages:
       * NumPy, SciPy, Matplotlib
       * NetworkX
       * MDTraj
   - Make scripts executable:
       chmod +x s1build s2gcmc s3md s4postanalysis

2. Build the Composite
   Run:
       ./s1build
   This constructs MOF/polymer composites using the MOF folders 
   (ALFFIVE, CALF20, MIL53, Zrfcufum) and polymers in all_pols/.

3. CO₂ Adsorption (GCMC)
   Run:
       ./s2gcmc
   Performs adsorption simulations with CADSS, generating 
   isotherms and adsorption states.

4. Molecular Dynamics
   Run:
       ./s3md
   Launches LAMMPS simulations of adsorbed CO₂ molecules.

5. Post-Analysis
   Run:
       ./s4postanalysis
   Performs:
     - Structural analysis (atomic density, RDFs, Delaunay tessellation)
     - Textural analysis (pore size distribution, void fraction, free-pore mapping)
     - Graph theory analysis (connectivity, betweenness, eccentricity)
     - Dynamics (diffusion coefficients, angular reorientation)
   Results are saved in the corresponding folders.

----------------------------------------------------------
Repository Structure
----------------------------------------------------------
├── s1build          # Build the composite system
├── ALFFIVE/         # Example MOF folder
├── CALF20/          # Example MOF folder
├── Zrfcufum/        # Example MOF folder
├── MIL53/           # Example MOF folder
├── all_pols/        # Polymer libraries
├── s2gcmc           # CO₂ adsorption with CADSS
├── s3md             # Molecular Dynamics with LAMMPS
├── inputs_msd/      # Input files for MSD/dynamics analysis
├── s4postanalysis   # Post-analysis (structure, texture, graph theory)
└── script/          # Utility scripts

----------------------------------------------------------
Quick Example
----------------------------------------------------------
# Build composite using MIL53 + polymer
./s1build MIL53

# Run CO₂ adsorption
./s2gcmc MIL53

# Run MD
./s3md MIL53

# Post-analysis
./s4postanalysis MIL53

This sequence builds, simulates, and analyzes one system end-to-end.

----------------------------------------------------------
Notes
----------------------------------------------------------
- Edit input parameters in inputs_msd/ before running the scripts.
- You can run each step independently, but the recommended order is s1 → s2 → s3 → s4.
"""
