function [connectome] = normalize_roisize_mrtrix(filename,roi_size,threshold)
% This function normalizes the connectivity matrix by the sum of volumes of
% all ROIs. As input, besides the connectivity matrix the list of ROI sizes
% in voxel count is given.

connectome=importdata(filename);

size_roi=importdata(roi_size);%.*8; %get size of rois in voxels and multiply by volume per voxel
mean_volume=mean(size_roi);

%Normalize MRTrix connectome with mean size of ROIs
n_nodes=length(connectome);
for n=1:n_nodes
    for m=1:n_nodes
        if n~=m
            connthreshold=threshold*2/(size_roi(n)+size_roi(m));
            if connectome(n,m)>connthreshold
                connectome(n,m)=connectome(n,m)*mean_volume/10000000;%Normalized by sizes of rois and number of streamlines
                %connectome(n,m)=connectome(n,m)/10000000/(2/(size_roi(n)+size_roi(m)));%Normalized by number of streamlines
                %connectome(n,m)=connectome(n,m)*mean_volume;
            else 
                connectome(n,m)=0;
            end
        end
    end
end
end

