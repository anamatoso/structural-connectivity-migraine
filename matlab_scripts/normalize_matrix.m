function [connectome] = normalize_matrix(filename,roi_size,threshold,algorithm,norm)
% This function normalizes the connectivity matrix by the sum of volumes of
% all ROIs. As input, besides the connectivity matrix the list (txt) of ROI sizes
% in voxel count is given.

connectome=load(filename);

size_roi=load(roi_size); %get size of rois in voxels 
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
                connthreshold=threshold*2/(size_roi(n)+size_roi(m)); %mrtrix comes with normalization by the volume of ROIs done
            end
            
            if connectome(n,m)>connthreshold
                if algorithm=='fsl'
                    if norm==1
                        connectome(n,m)=(connectome(n,m)/(5000*sum_volume))*2*mean_volume/(size_roi(n)+size_roi(m)); %Normalize by sizes of rois and number of streamlines
                    elseif norm==2
                        connectome(n,m)=(connectome(n,m))*2*mean_volume/(size_roi(n)+size_roi(m));%Normalize just by size of ROIs
                    elseif norm==3
                        connectome(n,m)=connectome(n,m)/(5000*sum_volume);%Normalize just by number of streamlines
                    end
                else
                    if norm==1
                        connectome(n,m)=connectome(n,m)*mean_volume/10000000;%Normalize by sizes of rois and number of streamlines
                    elseif norm==2
                        connectome(n,m)=connectome(n,m)*mean_volume;%Normalize just by size of ROIs
                    elseif norm==3
                        connectome(n,m)=connectome(n,m)/10000000/(2/(size_roi(n)+size_roi(m)));%Normalize by number of streamlines
                    else
                        connectome(n,m)=connectome(n,m)/(2/(size_roi(n)+size_roi(m)));%Not normalize
                    end
                end
            else 
                connectome(n,m)=0;
            end
        end
    end
end

% Make sure matrix is symmetrical (FSL)
if ~issymmetric(connectome)
    connectome=(connectome+connectome');
end

end

