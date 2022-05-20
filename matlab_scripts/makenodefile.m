function [matrix] = makenodefile(color_size)

coord=importdata('aal116_MNIcoord.txt');
labels=get_label_nodes("AAL116_labels.txt")';

matrix=[coord color_size labels];

end

