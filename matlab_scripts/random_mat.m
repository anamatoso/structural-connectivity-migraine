function [randmat] = random_mat(mat,nperm)
nnodes=length(mat);
randmat=mat;
for i=1:nperm
    %nodes randomly chosen
    n11=randi(nnodes);n12=randi(nnodes);
    n21=randi(nnodes);n22=randi(nnodes);
    
    while n11==n12 || n21==n22
        n11=randi(nnodes);n12=randi(nnodes);
        n21=randi(nnodes);n22=randi(nnodes);
    end
    val1=randmat(n11,n12);
    val2=randmat(n21,n22);
    
    randmat(n11,n12) = val2;
    randmat(n12,n11) = val2;
    
    randmat(n21,n22) = val1;
    randmat(n22,n21) = val1;
end
end

