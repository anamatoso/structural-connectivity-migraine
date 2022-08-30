function [connectome] = normalize_matrix(filename,roi_size,threshold,algorithm,norm)
% This function normalizes the connectivity matrix by the sum of volumes of
% all ROIs. As input, besides the connectivity matrix the list of ROI sizes
% in voxel count is given.

connectome=importdata(filename);

size_roi=importdata(roi_size);%.*8; %get size of rois in voxels and multiply by volume per voxel
mean_volume=mean(size_roi);
sum_volume=sum(size_roi);


%Normalize connectome
n_nodes=length(connectome);
for n=1:n_nodes
    for m=1:n_nodes
        if n~=m
            if algorithm=='fsl'
                connthreshold=threshold;
            else
                connthreshold=threshold*2/(size_roi(n)+size_roi(m));
            end
            
            if connectome(n,m)>connthreshold
                if algorithm=='fsl'
                    if norm==1
                        connectome(n,m)=(connectome(n,m)/(5000*sum_volume))*2*mean_volume/(size_roi(n)+size_roi(m)); %Normalized by sizes of rois and number of streamlines
                    else
                        connectome(n,m)=(connectome(n,m))*2*mean_volume/(size_roi(n)+size_roi(m));%Normalized just by size of ROIs
                    end
                    %connectome(n,m)=connectome(n,m)/(5000*sum_volume);%Normalized just by number of streamlines
                else
                    if norm==1
                        connectome(n,m)=connectome(n,m)*mean_volume/10000000;%Normalized by sizes of rois and number of streamlines
                    else
                        connectome(n,m)=connectome(n,m)*mean_volume;%Normalized just by size of ROIs
                    end
                    %connectome(n,m)=connectome(n,m)/10000000/(2/(size_roi(n)+size_roi(m)));%Normalized by number of streamlines
                end
            else 
                connectome(n,m)=0;
            end
        end
    end
end

% Make sure it is symmetrical
if ~issymmetric(connectome)
    connectome=(connectome+connectome')/2;
end

if norm==2
% Add normalization of the sum of fibers that were actually counted
    connectome=connectome./(sum(sum(connectome)));
end

end

