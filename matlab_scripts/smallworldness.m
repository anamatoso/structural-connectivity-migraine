function [S] = smallworldness(mat)
% This function calculates the small worldness of a given connectivity
% matrix. Since the clustering coefficient and the characteristic path
% length are already calculated, they are put as input as well.

len_mat=1./mat;                 % conection-length matrix
d_mat= distance_wei(len_mat);
C=mean(clustering_coef_wu(weight_conversion(mat, 'normalize'))); % clustering coefficient
[L,~]=charpath(d_mat);

% falta calcular equivalent nodes and edges de mat
%...


% nodes=116;
% edges=sum(sum(mat));
% maximum=116;minimum=1;
% s=[];t=[];
% for i=1:edges
%     s = [s minimum + (maximum-minimum).*rand(1,1)];
%     t = [t minimum + (maximum-minimum).*rand(1,1)];
% end
% G = graph(s, t, [], nodes);
% random_mat=full(adjacency(G));

% falta Go back to scale of mat
% ...

random_mat=randmio_und(mat,100);
len_mat_rand=1./random_mat;                 % conection-length matrix
d_mat_rand= distance_wei(len_mat_rand);
Crand=mean(clustering_coef_wu(weight_conversion(random_mat, 'normalize'))); % clustering coefficient
[Lrand,~]=charpath(d_mat_rand);

n=C/Crand;
d=L/Lrand;
S=n/d;
end

