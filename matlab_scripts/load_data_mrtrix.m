function [connectomes] = load_data_mrtrix(folder,pattern,folder_roi_sizes)
% This function load the data from a given folder which have a given
% pattern

F = folder;
filePattern = fullfile(F, pattern); % Change to whatever pattern you need.
theFiles = dir(filePattern);
baseFileName = theFiles(1).name;
fullFileName = fullfile(theFiles(1).folder, baseFileName);
c=importdata(fullFileName);
connectomes=zeros(length(c),length(c),length(theFiles));

for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    roi_file=strcat(folder_roi_sizes,"/",baseFileName(1:length(baseFileName)-34),"_roi_size_intersect.txt"); %17
    fullFileName = fullfile(theFiles(k).folder, baseFileName);
    c=normalize_roisize_mrtrix(fullFileName,roi_file);
    connectomes(:,:,k)=c;
end
end

