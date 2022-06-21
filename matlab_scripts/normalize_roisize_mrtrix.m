function [connectome] = normalize_roisize_mrtrix(filename,roi_size)
% This function normalizes the connectivity matrix by the sum of volumes of
% all ROIs. As input, besides the connectivity matrix the list of ROI sizes
% in voxel count is given.

connectome=importdata(filename);

size_roi=importdata(roi_size);%.*8; %get size of rois in voxels and multiply by volume per voxel
sum_volume=sum(size_roi);
%Normalize MRTrix connectome with sizes of ROIs
n_nodes=length(connectome);
for n=1:n_nodes
    for m=1:n_nodes
        if n~=m
            connectome(n,m)=connectome(n,m)*sum_volume;
        end
    end
end
end

