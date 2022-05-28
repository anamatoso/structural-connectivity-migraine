function [matrix] = makenodefile(nodecoord,labels,color_size)

coord=importdata(nodecoord);
labels=labels';

matrix=[coord color_size labels];

end

