function [new_matrix] = remap_matrix(matrix,idx_map)

[n_nodes,~]=size(matrix);
n_newnodes=length(unique(idx_map));
new_matrix=zeros(n_newnodes,n_newnodes);

for i=2:n_nodes-1
    for j=i+1:n_nodes
        new_i= idx_map(i);
        new_j= idx_map(j);
        if new_i~=new_j
            new_matrix(new_i,new_j)=new_matrix(new_i,new_j)+matrix(i,j);
            new_matrix(new_j,new_i)=new_matrix(new_j,new_i)+matrix(j,i);
        end
    end
end
end

