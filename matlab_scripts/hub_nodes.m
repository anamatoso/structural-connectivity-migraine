function [I] = hub_nodes(conmat)
len_mat=1./conmat;                 % conection-length matrix
BC=betweenness_wei(len_mat)'; 

P = prctile(BC,80);
I = find(BC >= P);

end

