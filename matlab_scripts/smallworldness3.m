function [S] = smallworldness3(mat)
% This function calculates the small worldness of a given connectivity
% matrix.
nnodes=length(mat);

len_mat=1./mat;                 % conection-length matrix
d_mat= distance_wei(len_mat);
C=mean(clustering_coef_wu(weight_conversion(mat, 'normalize'))); % clustering coefficient
[L,~]=charpath(d_mat,0,0);

% create random matrix
randmat = zeros(size(mat));
for i=1:nnodes-1
    for j=i+1:nnodes
        m=1;n=1;
        while (m==n) && randmat(m,n)~=0
            m=randi(nnodes);n=randi(nnodes);
        end
        randmat(m,n)=mat(i,j);
        randmat(n,m)=mat(i,j);
    end
end



len_mat_rand=1./randmat;                 % conection-length matrix
d_mat_rand= distance_wei(len_mat_rand);
Crand=mean(clustering_coef_wu(weight_conversion(randmat, 'normalize'))); % clustering coefficient
[Lrand,~]=charpath(d_mat_rand,0,0);

n=C/Crand;
d=L/Lrand;
S=n/d;
end

