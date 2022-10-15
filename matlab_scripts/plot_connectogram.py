import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
#pip3 install -U mne-connectivity
import mne
from mne.datasets import sample
from mne.minimum_norm import apply_inverse_epochs, read_inverse_operator
from mne.viz import circular_layout
from mne_connectivity import spectral_connectivity_epochs
from mne_connectivity.viz import plot_connectivity_circle
from tkinter.filedialog import askopenfilename


filename = askopenfilename() #input("What matrix file do you want to plot? ")
connectivity = sio.loadmat(filename)
connectivity=connectivity['matrix']

number_edges=int(np.count_nonzero(connectivity)/2)

#separate_hemispheres = input("Do you yant to separate both hemispheres?[0/1] ")
separate_hemispheres="1"

if separate_hemispheres=="1": 
    separate_hemispheres=True
else:
    separate_hemispheres=False
# Get labels of nodes
label_names = open("AAL116_labels.txt","r").read().splitlines()

# Create one group for left hemisphere and another for right hemisphere
lh_labels = [name for name in label_names if name.endswith('-L')]

label_ypos = list()
for name in lh_labels:
    idx = label_names.index(name)
    ypos = label_names[idx]
    label_ypos.append(ypos)

lh_labels = [label for (yp, label) in (zip(label_ypos, lh_labels))]

rh_labels = [name for name in label_names if name.endswith('-R')]

lh_labels=label_names[::2] 
rh_labels=label_names[1::2]

# Save the plot order and create a circular layout
node_order = list()
node_order.extend(lh_labels[::-1])  # reverse the order
node_order.extend(rh_labels)



if separate_hemispheres:
    group_boundaries_final=[13,len(label_names) / 2,len(label_names)-13]
    node_order_final=node_order
    #print(separate_hemispheres)
else:
    group_boundaries_final=None
    node_order_final=label_names
    print(separate_hemispheres)

node_angles = circular_layout(label_names, node_order_final, start_pos=90, group_boundaries=group_boundaries_final)

node_colors=[['blueviolet']*2,['magenta']*26,['red']*2,['lightcoral']*6,['red']*6,['darkorange']*4,['yellowgreen']*2,['darkorange']*6,['green']*2,['blueviolet']*2,['cyan']*10,['blueviolet']*2,['maroon']*6,['red']*2,['green']*12,['gray']*26]
node_colors = [item for sublist in node_colors for item in sublist]

fig, ax = plt.subplots(1,1,figsize=(20, 20), facecolor='white',subplot_kw=dict(projection="polar"))
plot_connectivity_circle(connectivity, label_names, n_lines=number_edges,node_colors=node_colors, node_angles=node_angles,facecolor='white', textcolor='black', node_edgecolor='white',ax=ax, colormap='bwr', colorbar_pos=(-0.1, 0.1),padding=3)
fig.tight_layout()
if os.path.isfile("connectome.png"):
    os.remove("connectome.png")
fig.savefig("connectome.png", facecolor=fig.get_facecolor())
