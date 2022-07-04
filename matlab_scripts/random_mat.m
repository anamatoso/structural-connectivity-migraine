function [randmat] = random_mat(mat,nperm)

randmat=mat;
for i=1:nperm
    %nodes randomly chosen
    n11=randi(116);n12=randi(116);
    n21=randi(116);n22=randi(116);
    
    while n11==n12 || n21==n22
        n11=randi(116);n12=randi(116);
        n21=randi(116);n22=randi(116);
    end
    val1=randmat(n11,n12);
    val2=randmat(n21,n22);
    
    randmat(n11,n12) = val2;
    randmat(n12,n11) = val2;
    
    randmat(n21,n22) = val1;
    randmat(n22,n21) = val1;
end
end

