function [] = plotconmat(filename,nstreamlines,t)
% This function plots the connectivity matrix by firstly normalizing it
% by the number of streamlines.

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

