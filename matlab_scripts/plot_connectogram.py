import sys
import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
# pip3 install -U mne-connectivity
from mne.viz import circular_layout
from mne_connectivity.viz import plot_connectivity_circle

# Example to plot matrix.mat (in which the variable in matlab is called matrix) separating the hemispheres (1). (
# Using AAL116) python3 plot_connectogram.py matrix.mat 1


# Get information from matrix file
filename = sys.argv[1]
connectivity = sio.loadmat(filename)
connectivity = connectivity['matrix']

# Get number of edges in the connectogram and whether the user want to separate hemispheres
number_edges = int(np.count_nonzero(connectivity) / 2)  # since the matrix is symmetric, there is twice the number of edges we will want in our graph
separate_hemispheres = int(sys.argv[2])  # second argument in function

# Get labels of nodes
label_names = open("AAL116_labels.txt", "r").read().splitlines()

# Define spaces between nodes if the user want to separate the hemispheres
if separate_hemispheres:

    # Create one group for left hemisphere and another for right hemisphere
    lh_labels = [name for name in label_names if name.endswith('-L')]

    label_ypos = list()
    for name in lh_labels:
        idx = label_names.index(name)
        ypos = label_names[idx]
        label_ypos.append(ypos)

    lh_labels = label_names[::2]
    rh_labels = label_names[1::2]

    # Save the plot order and create a circular layout
    node_order = list()
    node_order.extend(lh_labels[::-1])  # reverse the order
    node_order.extend(rh_labels)

    # Define variables for the circular layout
    group_boundaries_final = [13, len(label_names) / 2, len(label_names) - 13]
    node_order_final = node_order

else:
    # Define variables for the circular layout if the user does not want to separate hemispheres
    group_boundaries_final = None
    node_order_final = label_names

# Determine the nodes' angles in the connectogram
node_angles = circular_layout(label_names, node_order_final, start_pos=90, group_boundaries=group_boundaries_final)

# Assign colors depending on the region
node_colors = [['mediumpurple'] * 2, ['steelblue'] * 26, ['orangered'] * 2, ['salmon'] * 6, ['orangered'] * 6,
               ['darkorange'] * 4, ['darkorange'] * 2, ['darkorange'] * 6, ['yellowgreen'] * 2, ['mediumpurple'] * 2,
               ['orchid'] * 10, ['mediumpurple'] * 2, ['brown'] * 6, ['orangered'] * 2, ['yellowgreen'] * 12,
               ['gray'] * 26]
node_colors = [item for sublist in node_colors for item in sublist]

# Create plot
fig, ax = plt.subplots(1, 1, figsize=(20, 20), facecolor='white', subplot_kw=dict(projection="polar"))
plot_connectivity_circle(connectivity, label_names, n_lines=number_edges, node_colors=node_colors,
                         node_angles=node_angles, facecolor='white', textcolor='black', node_edgecolor='white', ax=ax,
                         colormap='bwr', colorbar_pos=(-0.1, 0.1), padding=3)
fig.tight_layout()
