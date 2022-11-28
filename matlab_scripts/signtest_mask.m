function [significance_mask] = signtest_mask(matrix)
% This function creates a significance mask where 1 means that the connection is significant and 0 otherwise given a sign-test.

[n_nodes,~,~]=size(matrix);
significance_mask=zeros(n_nodes,n_nodes);
n_comparisons=(n_nodes*n_nodes-n_nodes)/2;

for i=1:n_nodes
    for j=1:n_nodes
        if i~=j
            data=squeeze(matrix(i,j,:));
            k=kstest(data);
            if k==1
                [p,h] = signtest(data,0,'Tail','right','alpha',0.05/n_comparisons);
                
            else
                [h,p] = ttest2(data,0,'Tail','right','alpha',0.05/n_comparisons);
            end
            if isnan(h)
                h=0;
            end
            significance_mask(i,j)=h;
            
        end
    end
end
end


