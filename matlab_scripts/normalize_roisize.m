function [connectome] = normalize_roisize(filename,roi_size)


connectome=importdata(filename);

size_roi=importdata(roi_size).*8; %get size of rois in voxels and multiply by volume/voxel

%Normalize FSL connectome with sizes of ROIs
for n=1:116
    for m=1:116
        if n~=m
            connectome(n,m)=connectome(n,m)/(size_roi(n)*size_roi(m));
        end
        
    end
end
end

