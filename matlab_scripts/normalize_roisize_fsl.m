function [connectome] = normalize_roisize_fsl(filename,roi_size)
% This function normalizes the connectivity matrix by the volume of the
% ROIs.
% As input, besides the connectivity matrix the list of ROI sizes in voxel
% count is given.

connectome=importdata(filename);

size_roi=importdata(roi_size).*8; %get size of rois in voxels and multiply by volume per voxel

%Normalize FSL connectome with sizes of ROIs
for n=1:116
    for m=1:116
        if n~=m
            connectome(n,m)=connectome(n,m)/(size_roi(n)+size_roi(m));
        end
    end
end
end
