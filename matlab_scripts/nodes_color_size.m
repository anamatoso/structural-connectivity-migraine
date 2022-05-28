function [table] = nodes_color_size(pvalues,difference,threshold,atlas_labels)

% size: small (1) for insignifican nodes, and big (2) for significant nodes
% color: 1 for insignificant, 3 for positive, 2 for negative
n_nodes=length(atlas_labels);
table=zeros(n_nodes,2);

for i=1:n_nodes
    
    if pvalues(i)<threshold %"significant" diference
        table(i,2)=2;
        if difference(i)==-1
            table(i,1)=3;
        else
            table(i,1)=2;
        end
    else %not "significant" diference
        table(i,2)=1;
        table(i,1)=1;
    end
    
end
end