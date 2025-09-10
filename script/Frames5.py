#!/usr/bin/env python
# coding: utf-8

# Importing required libraries
import numpy as np
import mdtraj as md
import networkx as nx
import itertools
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import interp1d
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt
import time
import math
import warnings
from string import digits
from scipy.spatial.distance import pdist
from collections import defaultdict
import meshio
import pandas as pd
from scipy.spatial import cKDTree
import community
import sys

# Ignore warnings
warnings.filterwarnings("ignore")

# Function to remove digits from a string
def remove_digits(input_str):
    return input_str.translate(str.maketrans('', '', digits))

# Read atomic radii from files and write to new files
mass = np.loadtxt("/home/adiaz/scripts/mass", dtype='str', delimiter=" ")
rad = np.loadtxt("/home/adiaz/scripts/radius", dtype='str', delimiter=" ")
lm = np.loadtxt("../listm.txt", dtype='str')
lp = np.loadtxt("../listp.txt", dtype='str')

with open('mass.mass', 'w+') as f_new1, open('rad.rad', 'w+') as f_new2:
    for val in lm:
        val = remove_digits(val)[:2].replace("_", "")
        for j in range(len(mass)):
            if val == mass[j, 1]:
                f_new1.write("%s\t%s\n" % (mass[j, 1], mass[j, 3]))
                num = float(rad[j, 3]) / 100
                f_new2.write("%s\t%f\n" % (mass[j, 1], num))
    for val in lp:
        val = remove_digits(val)[:1].replace("_", "").capitalize()
        for j in range(len(mass)):
            if val == mass[j, 1]:
                f_new1.write("%s\t%s\n" % (mass[j, 1], mass[j, 3]))
                num = float(rad[j, 3]) / 100
                f_new2.write("%s\t%f\n" % (mass[j, 1], num))

# Load trajectory and topology
traj1 = md.load('../wrapped.dcd', top='../first-final.pdb')  # Trajectory read in nanometer unit
traj1.xyz = 10 * traj1.xyz  # Change the unit back to angstrom
topology = traj1.topology  # Topology information of the system
no_atoms = traj1.n_atoms  # Total number of atoms
no_frames = 11 # Total number of frames
space = 1 # Total number of frames
box_size = 10*traj1.unitcell_lengths[0]
print (box_size)
if box_size[2]>= 1000:
    sys.exit("The box size in the z-axis is greater than or equal to 1000. Exiting the script.")    


# Set parameters for grid generation
size = 0.6 # 0.8
g = 0.8 # 0.8
g2 = g * 3.0
mind = 2 * size
maxsize = 10
gs=g
gs3=gs*2
sizes=size
with open("num_dots.txt", "w") as f:
    f.write(str(int(box_size[0]*box_size[1]/gs)))

# Initialize a dictionary to store results for each frame and label
results_per_frame_system1 = defaultdict(list)
results_per_frame_system2 = defaultdict(list)
results_per_frame = defaultdict(list)
# Iterate through frames
for frame in range(0, no_frames, space):
    # Get the volume of the unit cell for the current frame
    vol = traj1.unitcell_volumes[frame] * 1000  # Volume in angstrom^n
    molecule_atoms = topology.select(f"resSeq {444}")
    com = md.compute_center_of_mass(traj1[frame].atom_slice(molecule_atoms)).flatten()
    translation = np.array([0, 0, (box_size[2] / 2 - com[2])])
    traj1.xyz[frame] += translation
    traj1.xyz[frame] -= np.floor(traj1.xyz[frame] / box_size) * box_size
    lowest_z_atom = traj1.atom_slice(molecule_atoms).xyz[frame, :, 2].argmin()
    highest_z_atom = traj1.atom_slice(molecule_atoms).xyz[frame, :, 2].argmax()
    lowest_z_coordinate = traj1.xyz[frame, lowest_z_atom, 2]
    highest_z_coordinate = traj1.xyz[frame, highest_z_atom, 2]

    zmin = round(lowest_z_coordinate) - 30
    zmax = round(highest_z_coordinate) + 30

    # Read atomic radii from file
    radii = {}
    with open("rad.rad", "r") as f:
        for line in f:
            atom_type, radius = line.strip().split()
            radii[atom_type] = float(radius)

    # Get coordinates from trajectory
    coords = []
    for atom in topology.atoms:
        atom_type = atom.element.symbol
        coord = traj1.xyz[frame][atom.index]
        radius = radii.get(atom_type, 1.0)
        coords.append((coord, radius))

    # Filter coordinates and radii based on zmin and zmax
    coord_array = []
    radius_array = []
    for coord, radius in coords:
        if zmin <= coord[2] <= zmax:
            coord_array.append(coord)
            radius_array.append(radius)

    arr = np.array(coord_array)
    radius_array = np.array(radius_array)

    # Define parameters
    num_bins = 100  # Number of bins along the z-axis
    extra = 2
    bin_counts = np.zeros(num_bins)
    z_coordinates = arr[:, 2]  # Extract z-coordinates from the filtered coordinates
    counts, bin_edges = np.histogram(z_coordinates, bins=num_bins, range=(zmin, zmax))
    bin_counts += counts
    bin_width = (zmax - zmin) / num_bins
    density = bin_counts / (bin_width * np.mean(box_size))
    z_density = [count / bin_width for count in bin_counts]
    z_values = np.linspace(zmin, zmax, num_bins)
    normalized_density = np.array(z_density) / sum(z_density)
    results_per_frame[frame].append(density)
    increased_density=density

    # Find the index corresponding to middle_z_coordinate
    middle_z_coordinate = (lowest_z_coordinate + highest_z_coordinate) / 2
    middle_z_index = np.abs(z_values - middle_z_coordinate).argmin()
    lowest_z_index = np.abs(z_values - lowest_z_coordinate).argmin()
    highest_z_index = np.abs(z_values - highest_z_coordinate).argmin()

    # Split the density into left and right halves
    left_density = increased_density[:middle_z_index]
    right_density = increased_density[middle_z_index:]

    # Find the indices of the local minima on the left side of the left maximum
    left_minima_indices = argrelextrema(left_density, np.less)[0]
    left_min_density_z = z_values[:middle_z_index][left_minima_indices]
    left_min_density = left_density[left_minima_indices]
    left_min_distances_left_max = (left_min_density_z - lowest_z_coordinate)
    left_min_distances_left_max[left_min_density_z >= lowest_z_coordinate] = np.inf
    left_min_index = np.argmin(np.abs(left_min_distances_left_max))
    left_min_density_value = left_min_density[left_min_index]
    left_min_density_z_value = left_min_density_z[left_min_index]


    # Find the indices of the local minima on the right side of the left maximum
    right_minima_indices = argrelextrema(right_density, np.less)[0]
    right_min_density_z = z_values[middle_z_index:][right_minima_indices]
    right_min_density = right_density[right_minima_indices]
    right_min_distances_left_max = (right_min_density_z - highest_z_coordinate)
    right_min_distances_left_max[right_min_density_z <= highest_z_coordinate] = np.inf
    right_min_index = np.argmin(np.abs(right_min_distances_left_max))
    right_min_density_value = right_min_density[right_min_index]
    right_min_density_z_value = right_min_density_z[right_min_index]

    # Find the indices of the local maxima on the left side of the left minimum
    left_maxima_indices = argrelextrema(left_density, np.greater)[0]
    left_max_density_z = z_values[:middle_z_index][left_maxima_indices]
    # Find the index of the closest maximum to the right of the left minimum
    left_closest_max_right_min_distance = left_max_density_z - left_min_density_z_value
    left_closest_max_right_min_distance[left_closest_max_right_min_distance <= 0] = np.inf
    left_closest_max_right_min_index = np.argmin(left_closest_max_right_min_distance)
    left_closest_max_right_min_density_value = left_density[left_maxima_indices[left_closest_max_right_min_index]]
    left_closest_max_right_min_density_z_value = left_max_density_z[left_closest_max_right_min_index]+7
    # Find the index of the closest maximum to the left of the left minimum
    left_closest_max_left_min_distance = left_min_density_z_value - left_max_density_z
    left_closest_max_left_min_distance[left_closest_max_left_min_distance <= 0] = np.inf
    left_closest_max_left_min_index = np.argmin(left_closest_max_left_min_distance)
    left_closest_max_left_min_density_value = left_density[left_maxima_indices[left_closest_max_left_min_index]]
    left_closest_max_left_min_density_z_value = left_max_density_z[left_closest_max_right_min_index]-25

    # Find the indices of the local maxima on the right side of the right minimum
    right_maxima_indices = argrelextrema(right_density, np.greater)[0]
    right_max_density_z = z_values[middle_z_index:][right_maxima_indices]
    # Find the index of the closest maximum to the right of the left minimum
    right_closest_max_left_min_distance = right_min_density_z_value - right_max_density_z
    right_closest_max_left_min_distance[right_closest_max_left_min_distance <= 0] = np.inf
    right_closest_max_left_min_index = np.argmin(right_closest_max_left_min_distance)
    right_closest_max_left_min_density_value = right_density[right_maxima_indices[right_closest_max_left_min_index]]
    right_closest_max_left_min_density_z_value = right_max_density_z[right_closest_max_left_min_index]-7
    # Find the index of the closest maximum to the right of the left minimum
    right_closest_max_right_min_distance = right_max_density_z - right_min_density_z_value
    right_closest_max_right_min_distance[right_closest_max_right_min_distance <= 0] = np.inf
    right_closest_max_right_min_index = np.argmin(right_closest_max_right_min_distance)
    right_closest_max_right_min_density_value = right_density[right_maxima_indices[right_closest_max_right_min_index]]
    right_closest_max_right_min_density_z_value = right_max_density_z[right_closest_max_left_min_index]+25
    # Save grid to file
    with open(f"reference_right{frame}.xyz", "w") as f:
        f.write(str(highest_z_coordinate))
    with open(f"reference_left{frame}.xyz", "w") as f:
        f.write(str(lowest_z_coordinate))

    dots = []
    for x, y, z in itertools.product(
        np.arange(0, box_size[0], gs),
        np.arange(0, box_size[1], gs),
        np.arange(left_closest_max_left_min_density_z_value, left_closest_max_right_min_density_z_value, gs)
    ):
        dist = np.sqrt(np.sum((arr[:, :3] - np.array([x, y, z])) ** 2, axis=1))
        in_range = np.logical_and(dist <= radius_array + sizes, np.logical_or(np.abs(x - arr[:, 0]) <= mind, np.abs(y - arr[:, 1]) <= mind, np.abs(z - arr[:, 2]) <= mind))
        if not np.any(in_range):
            dots.append(np.array([x, y, z]))

    for x, y, z in itertools.product(
        np.arange(0, box_size[0], gs),
        np.arange(0, box_size[1], gs),
        np.arange(right_closest_max_left_min_density_z_value, right_closest_max_right_min_density_z_value, gs)
    ):
        dist = np.sqrt(np.sum((arr[:, :3] - np.array([x, y, z])) ** 2, axis=1))
        in_range = np.logical_and(dist <= radius_array + sizes, np.logical_or(np.abs(x - arr[:, 0]) <= mind, np.abs(y - arr[:, 1]) <= mind, np.abs(z - arr[:, 2]) <= mind))
        if not np.any(in_range):
            dots.append(np.array([x, y, z]))

    dots = np.unique(dots, axis=0)

    # Save grid to file
    with open(f"interphase_{frame}.xyz", "w") as f:
        for dot in dots:
            f.write(" ".join(str(x) for x in dot) + "\n")

    # Create two separate graphs based on z-coordinate
    G_filtered_z_above = nx.Graph()
    G_filtered_z_below = nx.Graph()
    
    # Load the VPSDPTS file into a graph
    pos = {}
    nodes_to_reverse = []  # Create a list to store nodes to be reversed
    xyz_filename = f"interphase_{frame}.xyz"
    
    with open(xyz_filename) as f:
        nodes = [tuple(map(lambda coord: round(float(coord), 2), line.split())) for line in f]
    
    for node in nodes:
        x, y, z = node
        pos[node] = (x, y, z)
    
        if z > middle_z_coordinate:
            G_filtered_z_above.add_node(node)
        else:
            G_filtered_z_below.add_node(node)
            #nodes_to_reverse.append(node)
    
    # Reverse the orientation of G_filtered_z_below in the z-direction
    #G_filtered_z_below.add_nodes_from((x, y, -z) for x, y, z in nodes_to_reverse)
    #pos.update((node, (x, y, -z)) for node, (x, y, z) in pos.items())
    
    # Save the graphs G_filtered_z_above and G_filtered_z_below to GraphML files
    nx.write_graphml(G_filtered_z_above, f"G_filtered_z_right{frame}.graphml")
    nx.write_graphml(G_filtered_z_below, f"G_filtered_z_left{frame}.graphml")
