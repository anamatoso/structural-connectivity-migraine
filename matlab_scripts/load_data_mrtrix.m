function [connectomes] = load_data_mrtrix(folder,pattern)
% This function load the data from a given folder which have a given
% pattern

F = folder;
filePattern = fullfile(F, pattern); % Change to whatever pattern you need.
theFiles = dir(filePattern);
connectomes=zeros(116,116,length(theFiles));
for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(theFiles(k).folder, baseFileName);
    c=importdata(fullFileName);
    connectomes(:,:,k)=c;
end
end

