import pandas as pd
import matplotlib.pyplot as plt
import os
import itertools as tee
import networkx as nx
import numpy as np
from itertools import combinations
from networkx.algorithms import community
from networkx.algorithms.similarity import graph_edit_distance
import seaborn as sns
from itertools import combinations, product
from scipy.spatial import cKDTree
from collections import Counter
import psutil
from matplotlib.cm import get_cmap
from scipy.interpolate import interp1d
from math import sqrt

fin = open('MonteCarloFinConfig.car', 'r')

frames = []
n_atoms = 0
n_frames = 1
frame_atoms = []
CO2_atoms = []
cell = []

rl = fin.readline
# line 1-4 is comment/title
dump=rl()
dump=rl()
dump=rl()
dump=rl()
#lines 5 contains PBC information
line = rl().strip()
cols=line.split()
if (cols[0] != 'PBC'):
    print ('check the file, line 5 should have PBC')
cell_a = float(cols[1])
cell_b = float(cols[2])
cell_c = float(cols[3])
cell_alpha = float(cols[4])
cell_beta = float(cols[5])
cell_gamma = float(cols[6])
cell = (cell_a,cell_b,cell_c,cell_alpha,cell_beta,cell_gamma)
print (cell)

DONE = False
reading_atom_sites=False

atomNr = 0
CO2Nr = 0
resNr = 1
resName="MOF"

while not DONE:
    line = rl().strip()
    cols=line.split()
    if (len(cols) == 9):
        atomNr += 1
        atomLabel,x,y,z,resNamex,resNrx,mmType,elType,charge = cols

        if (mmType == "O_e"):
            resNr += 1
        if (mmType == "CO2_C"):
            resNr += 1        
        if (mmType == "O_e" or mmType == "H_e"):
            resName = "HOH"
        elif (mmType == "CO2_C" or mmType == "CO2_O"):
            resName = "CO2"
            CO2Nr +=1
        elif (mmType == "N2_N"):
            resName = "N2"
        else:
            resName = "MOF"
        frame_atoms.append((atomLabel, float(x), float(y), float(z), resName, int(resNr), mmType, elType, float(charge)))
        if resName != "MOF":
            CO2_atoms.append((atomLabel, float(x), float(y), float(z), resName, int(resNr), mmType, elType, float(charge)))
    # if empty line, we're done
    if not line.strip():
        DONE = True
        #print ("DONE reading ", atomNr, "atoms" )
        fin.close()
        #otherwise just continue with the loop
        frames.append((cell, frame_atoms))
        #print("C- ", frames)


# Open the LAMMPS data file
with open("../final-interface.lmps", "r") as file:
    # Read the lines of the file
    lines = file.readlines()
    # Iterate over each line
    for line in lines:
        # Split the line into columns
        cols = line.split()

        # Check if the line contains "atoms"
        if "atoms" in line:
            Nratoms = int(cols[0])
        elif "atom types" in line:
            atom_types = int(cols[0])
        elif "bonds" in line:
            bonds = int(cols[0])
        elif "bond types" in line:
            bond_types = int(cols[0])
        elif "angles" in line:
            angles = int(cols[0])
        elif "angle types" in line:
            angle_types = int(cols[0])
        elif "dihedrals" in line:
            dihedrals = int(cols[0])
        elif "dihedral types" in line:
            dihedral_types = int(cols[0])
        elif "impropers" in line:
            impropers = int(cols[0])
        elif "improper types" in line:
            improper_types = int(cols[0])
        elif "xhi" in line:
            xlow = float(cols[0])
            xhig = float(cols[1])
        elif "yhi" in line:
            ylow = float(cols[0])
            yhig = float(cols[1])
        elif "zhi" in line:
            zlow = float(cols[0])
            zhig = float(cols[1])
    # Print the extracted information
    print("Number of Atoms:", Nratoms)
    print("Number of Atom Types:", atom_types)
    print("Number of Bonds:", bonds)
    print("Number of Bond Types:", bond_types)
    print("Number of Angles:", angles)
    print("Number of Angle Types:", angle_types)
    print("Number of Dihedrals:", dihedrals)
    print("Number of Dihedral Types:", dihedral_types)
    print("Number of Impropers:", impropers)
    print("Number of Improper Types:", improper_types)
    # List to store the lines before "Bonds"
    atoms_lines = []
    bonds_lines = []
    angles_lines = []
    dihedrals_lines = []
    impropers_lines = []
    
    # Check atom lines
    start_storing = False
    check = False
    for line in lines:
        if "Bonds" in line:
            break
        if start_storing:
            if check:
                if line.strip() == "":
                    break
                atoms_lines.append(line.strip().split())
            else:
                check =True
        if "Atoms" in line:
            start_storing = True
    sorted_atoms_lines = sorted(atoms_lines, key=lambda x: int(x[0]))

    # Check bond lines
    start_storing = False
    check = False
    for line in lines:
        if "Angles" in line:
            break
        if start_storing:
            if check:
                if line.strip() == "":
                    break
                bonds_lines.append(line.strip().split())
            else:
                check =True
        if "Bonds" in line:
            start_storing = True
    sorted_bonds_lines = sorted(bonds_lines, key=lambda x: int(x[0]))

    # Check angles lines
    start_storing = False
    check = False
    for line in lines:
        if "Dihedrals" in line:
            break
        if start_storing:
            if check:
                if line.strip() == "":
                    break
                angles_lines.append(line.strip().split())
            else:
                check =True
        if "Angles" in line:
            start_storing = True
    sorted_angles_lines = sorted(angles_lines, key=lambda x: int(x[0]))

    # Check dihedrals lines
    start_storing = False
    check = False
    for line in lines:
        if "Impropers" in line:
            break
        if start_storing:
            if check:
                if line.strip() == "":
                    break
                dihedrals_lines.append(line.strip().split())
            else:
                check =True
        if "Dihedrals" in line:
            start_storing = True
    sorted_dihedrals_lines = sorted(dihedrals_lines, key=lambda x: int(x[0]))

    # Check Impropers lines
    start_storing = False
    check = False
    for line in lines:
        if start_storing:
            if check:
                if line.strip() == "":
                    break
                impropers_lines.append(line.strip().split())
            else:
                check =True
        if "Impropers" in line:
            start_storing = True
    sorted_impropers_lines = sorted(impropers_lines, key=lambda x: int(x[0]))

    # Check Impropers lines
    start_storing = False
    check = False
    lowest_z = float('inf')
    highest_z = float('-inf')
    for line in lines:
        if "Bonds" in line:
            break
        if start_storing:
            if check:
                if line.strip() == "":
                    break
                columns = line.strip().split()
                # Check if the atom belongs to the 444 group (assuming the third column indicates the group)
                if columns[1] == '444':
                    z_coord = float(columns[6])
                    lowest_z = min(lowest_z, z_coord)
                    highest_z = max(highest_z, z_coord)
            else:
                check =True
        if "Atoms" in line:
            start_storing = True

    print (abs(lowest_z-zlow), abs(highest_z-zhig))
    if abs(lowest_z-zlow)<1 or abs(highest_z-zhig)<1:
        start_storing = False
        check = False
        lowest_z = float('inf')
        highest_z = float('-inf')
        print ("MOF is crossing the boundary")
        zmid=(zhig-zlow)/2+zlow
        for line in lines:
            if "Bonds" in line:
                break
            if start_storing:
                if check:
                    if line.strip() == "":
                        break
                    columns = line.strip().split()
                    # Check if the atom belongs to the 444 group (assuming the third column indicates the group)
                    if columns[1] == '444':
                        z_coord = float(columns[6])
                        if z_coord >zmid:
                            lowest_z = min(lowest_z, z_coord)
                        else:
                            highest_z = max(highest_z, z_coord)
                else:
                    check =True
            if "Atoms" in line:
                start_storing = True

    # Print the lowest and highest z-coordinate values
    print("Lowest z-coordinate:", lowest_z)
    print("Highest z-coordinate:", highest_z)

fout = open('CO2_all_dummies.lmps', 'w')
w = fout.write
w('LAMMPS data file via alejandro\n\n')

w(f'{Nratoms+CO2Nr} atoms\n')
w(f'{atom_types+4} atom types\n')
w(f'{bonds+int(CO2Nr/3*2)} bonds\n')
w(f'{bond_types+2} bond types\n')
w(f'{angles+int(CO2Nr/3)} angles\n')
w(f'{angle_types+2} angle types\n')
w(f'{dihedrals} dihedrals\n')
w(f'{dihedral_types} dihedral types\n')
w(f'{impropers} impropers\n')
w(f'{improper_types} improper types\n\n')

icell_a=0.0
icell_b=0.0
icell_c=0.0
adiff=cell_a-icell_a
bdiff=cell_b-icell_b
cdiff=cell_c-icell_c

w(f'{xlow} {xhig} xlo xhi\n')
w(f'{ylow} {yhig} ylo yhi\n')
w(f'{zlow} {zhig} zlo zhi\n\n')
xdiff=xhig-xlow
ydiff=yhig-ylow
zdiff=zhig-zlow

w('Atoms # full\n\n')
sgroup = max(int(line[1]) for line in sorted_atoms_lines)
satom_type = max(int(line[2]) for line in sorted_atoms_lines)

def move_atoms_inside_bounds(j,atoms, icell_a, cell_a, adiff, icell_b, cell_b, bdiff, icell_c, cell_c, cdiff, sgroup, satom_typei,O_pass):
    atom = list(atoms)
    bx, by, bz = 0, 0, 0  # Initialize box movement flags

    # Adjust atom position in x-direction
    if atom[1] < icell_a:
        bx = -1
        atom[1] += adiff
    elif atom[1] > cell_a:
        bx = 1
        atom[1] -= adiff

    # Adjust atom position in y-direction
    if atom[2] < icell_b:
        by = -1
        atom[2] += bdiff
    elif atom[2] > cell_b:
        by = 1
        atom[2] -= bdiff

    ## Adjust atom position in z-direction
    if atom[3] < icell_c:
        bz = -1
        atom[3] += cdiff
    elif atom[3] > cell_c:
        bz = 1
        atom[3] -= cdiff

    # Update atom tuple
    atom = tuple(atom)
    # add limits:
    lll=lowest_z-10#43.903500#+14#+icell_c/4
    llh=lowest_z-00#52.403500#+14#+icell_c/4
    hll=highest_z+00#90.403500#+14#+icell_c/4
    hlh=highest_z+10#98.903500#+14#+icell_c/4
    
    if lll<atom[3]<llh or hll<atom[3]<hlh:
        if atom[6] == 'CO2_C':
            moved_atoms=(f'{j} {atom[5]+sgroup-1} {satom_type+1} {atom[8]} {atom[1]} {atom[2]} {atom[3]} {bx} {by} {bz}\n')
            O_pass=2
    else:
        if atom[6] == 'CO2_C':
            moved_atoms=(f'{j} {atom[5]+sgroup-1} {satom_type+3} {atom[8]} {atom[1]} {atom[2]} {atom[3]} {bx} {by} {bz}\n')
    if atom[6] == 'CO2_O':
        if O_pass<=0: 
            moved_atoms=(f'{j} {atom[5]+sgroup-1} {satom_type+4} {atom[8]} {atom[1]} {atom[2]} {atom[3]} {bx} {by} {bz}\n')
        else:
            moved_atoms=(f'{j} {atom[5]+sgroup-1} {satom_type+2} {atom[8]} {atom[1]} {atom[2]} {atom[3]} {bx} {by} {bz}\n')
            O_pass-=1

#    if lll<atom[3]<llh or hll<atom[3]<hlh:
#        # Write the atom details in LAMMPS format
#        if atom[6] == 'CO2_C':
#            moved_atoms=(f'{j} {atom[5]+sgroup-1} {satom_type+1} {atom[8]} {atom[1]} {atom[2]} {atom[3]} {bx} {by} {bz}\n')
#        elif atom[6] == 'CO2_O':
#            moved_atoms=(f'{j} {atom[5]+sgroup-1} {satom_type+2} {atom[8]} {atom[1]} {atom[2]} {atom[3]} {bx} {by} {bz}\n')
#    else:
#        #print ('outside')
#        if atom[6] == 'CO2_C':
#            moved_atoms=(f'{j} {atom[5]+sgroup-1} {satom_type+3} {atom[8]} {atom[1]} {atom[2]} {atom[3]} {bx} {by} {bz}\n')
#        elif atom[6] == 'CO2_O':
#            moved_atoms=(f'{j} {atom[5]+sgroup-1} {satom_type+4} {atom[8]} {atom[1]} {atom[2]} {atom[3]} {bx} {by} {bz}\n')

    return (moved_atoms,O_pass)

O_pass = 0
for i, (atoms) in enumerate(frame_atoms):
    j=i+1
    atoms = tuple(atoms)
    if atoms[6] == 'CO2_C' or atoms[6] == 'CO2_O':
        # Call the function with updated O_pass
        moved_atoms_tuple, O_pass = move_atoms_inside_bounds(j, atoms, xlow, xhig, xdiff, ylow, yhig, ydiff, zlow, zhig, zdiff, sgroup, satom_type, O_pass)
        #moved_atoms = moved_atoms_tuple[0]  # Extract the moved_atoms string
        #print (moved_atoms_tuple,O_pass)
        w(moved_atoms_tuple)
        #w(move_atoms_inside_bounds(j, atoms, xlow, xhig, xdiff, ylow, yhig, ydiff, zlow, zhig, zdiff, sgroup, satom_type,O_pass))
    else:
        w(f'{j} {sorted_atoms_lines[i][1]} {sorted_atoms_lines[i][2]} {sorted_atoms_lines[i][3]} {sorted_atoms_lines[i][4]} {sorted_atoms_lines[i][5]} {sorted_atoms_lines[i][6]} {sorted_atoms_lines[i][7]} {sorted_atoms_lines[i][8]} {sorted_atoms_lines[i][9]}\n')

w('\n')
w('Bonds\n\n')
for atoms in sorted_bonds_lines:
    w(f'{atoms[0]} {atoms[1]} {atoms[2]} {atoms[3]}\n')
k=bonds
for i, (atoms) in enumerate(CO2_atoms):
    j=Nratoms+i+1
    if atoms[6] == 'CO2_C':
        k+=1
        w(f'{k} {bond_types+1} {j} {j+1}\n')
        k+=1
        w(f'{k} {bond_types+1} {j} {j+2}\n')
w('\n')
w('Angles\n\n')
for atoms in sorted_angles_lines:
    w(f'{atoms[0]} {atoms[1]} {atoms[2]} {atoms[3]} {atoms[4]}\n')
k=angles
for i, (atoms) in enumerate(CO2_atoms):
    j=Nratoms+i+1
    if atoms[6] == 'CO2_C':
        k+=1
        w(f'{k} {angle_types+1} {j+1} {j} {j+2}\n')
w('\n')
w('Dihedrals\n\n')
for atoms in sorted_dihedrals_lines:
    w(f'{atoms[0]} {atoms[1]} {atoms[2]} {atoms[3]} {atoms[4]} {atoms[5]}\n')
w('\n')
w('Impropers\n\n')
for atoms in sorted_impropers_lines:
    w(f'{atoms[0]} {atoms[1]} {atoms[2]} {atoms[3]} {atoms[4]} {atoms[5]}\n')
