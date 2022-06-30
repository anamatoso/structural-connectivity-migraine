function [connectomes] = load_data_mrtrix_OG(folder,pattern)
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
    fullFileName = fullfile(theFiles(k).folder, baseFileName);
    connectomes(:,:,k)=importdata(fullFileName);
end
end
