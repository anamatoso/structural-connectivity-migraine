function [randmat] = random_matrix(mat)
% create random matrix
randmat = zeros(size(mat));
nnodes=length(mat);
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
end