function [table] = nodes_color_size(pvalues,difference,threshold)
% This function creates a matrix that will function as input in BrainNet

% Size: small (=1) for insignifican nodes, and big (=3) for significant nodes
% Colour: 1 for insignificant, 2 for positive, 0 for negative

n_nodes=length(pvalues);
table=zeros(n_nodes,2); % color | size

for i=1:n_nodes
    if pvalues(i)<=threshold  % significant diference
        table(i,2)=3;
        if difference(i)==1
            table(i,1)=2;
        else
            table(i,1)=0;
        end
        
    else %not significant diference
        table(i,2)=1;
        table(i,1)=1;
    end
end
end