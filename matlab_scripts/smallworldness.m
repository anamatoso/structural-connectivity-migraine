function [S] = smallworldness(mat)
% This function calculates the small worldness of a given connectivity
% matrix. Since the clustering coefficient and the characteristic path
% length are already calculated, they are put as input as well.
mat=round(mat.*1e30);


len_mat=1./mat;                 % conection-length matrix
d_mat= distance_wei(len_mat);
C=mean(clustering_coef_wu(weight_conversion(mat, 'normalize'))); % clustering coefficient
[L,~]=charpath(d_mat,0,0);

random_mat=randmio_und(mat,100);
len_mat_rand=1./random_mat;                 % conection-length matrix
d_mat_rand= distance_wei(len_mat_rand);
Crand=mean(clustering_coef_wu(weight_conversion(random_mat, 'normalize'))); % clustering coefficient
[Lrand,~]=charpath(d_mat_rand,0,0);

n=C/Crand;
d=L/Lrand;
S=n/d;
end

