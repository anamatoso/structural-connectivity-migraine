function [significance_mask] = signtest_mask(a)
% a is matrix of size 116x116xn_people. Comes from load data.
% significance_mask is 1 is connection is significant and 0 otherwise
[n_nodes,~,~]=size(a);
significance_mask=zeros(n_nodes,n_nodes);
n_comparisons=(n_nodes*n_nodes-n_nodes)/2;

for i=1:n_nodes
    for j=1:n_nodes
        if i~=j
            data=squeeze(a(i,j,:));
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


