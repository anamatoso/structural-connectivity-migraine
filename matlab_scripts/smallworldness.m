function [S] = smallworldness(C,L,mat)
% This function calculates the small worldness of a given connectivity
% matrix. Since the clustering coefficient and the characteristic path
% length are already calculated, they are put as input as well.


% falta calcular equivalent nodes and edges de mat
%...


s = randi(nodes, edges, 1);
t = randi(nodes, edges, 1);
G = graph(s, t, [], nodes);
random_mat=full(adjacency(G));

% falta Go back to scale of mat
% ...

d_mat= distance_wei(random_mat);
[Lrand,~]=charpath(d_mat);
Ci=clustering_coef_wu(weight_conversion(random_mat, 'normalize'));
Crand=mean(Ci);

n=C/Crand;
d=L/Lrand;
S=n/d;
end

