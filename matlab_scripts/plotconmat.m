function [] = plotconmat(filename,nstreamlines,t)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if isfile("matrices data/"+filename)
    connectome=importdata(filename)./nstreamlines;
else
    connectome=filename./nstreamlines;
end

%connectome=importdata(filename);
imagesc(connectome)
title(t)
colormap jet
colorbar
end

