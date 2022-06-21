function [connectomes] = load_data_fsl(folder_matrix,pattern,folder_roi_sizes)
% This function load the data from a given folder which have a given
% pattern

F = folder_matrix;
filePattern = fullfile(F, pattern);
theFiles = dir(filePattern);
connectomes=zeros(116,116,length(theFiles));
for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    roi_file=strcat(folder_roi_sizes,"/",baseFileName(1:length(baseFileName)-46),"_roi_size.txt"); %17
    fullFileName = fullfile(theFiles(k).folder, baseFileName);
    c=normalize_roisize_fsl(fullFileName,roi_file);
    connectomes(:,:,k)=c;
end
end

