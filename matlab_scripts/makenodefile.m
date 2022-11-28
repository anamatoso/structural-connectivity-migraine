function [matrix] = makenodefile(nodecoord,labels,color_size)
% This functions creades th node file to be used in BrainNet

coord=importdata(nodecoord);
labels=labels';

matrix=[coord color_size labels];

end

