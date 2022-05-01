function [connectomes] = load_data(folder,pattern)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
F = folder;
filePattern = fullfile(F, pattern); % Change to whatever pattern you need.
theFiles = dir(filePattern);
connectomes=zeros(116,116,length(theFiles));
for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(theFiles(k).folder, baseFileName);
    % Now do whatever you want with this file name,
    % such as reading it in as an image array with imread()
    c=importdata(fullFileName);
    connectomes(:,:,k)=c;
end
end

