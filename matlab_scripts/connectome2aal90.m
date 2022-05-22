function [c] = connectome2aal90(connectomes)
c=cell(size(connectomes));
for i=1:length(connectomes)
    matrices=connectomes{i};
    matrices=matrices(1:90,1:90,:);
    c{i}=matrices;
end
end