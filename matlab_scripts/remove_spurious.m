function [significance_mask] = remove_spurious(group_matrix,mode,threshold)
% a is matrix of size 116x116xn_people. Comes from load data.
% significance_mask is 1 is connection is significant and 0 otherwise
% mode=1 is numer of fibers, mode=2 is percentil

[n_nodes,~,n_people]=size(group_matrix);
significance_mask=zeros(n_nodes,n_nodes);

if mode==1
    for p=1:n_people
        matrix=squeeze(group_matrix(:,:,p));
        for i=1:n_nodes-1
            for j=i:n_nodes
                if matrix(i,j)>threshold
                    significance_mask(i,j)=1;
                    significance_mask(j,i)=1;
                end
            end
        end
    end

else
    for p=1:n_people
        matrix=squeeze(group_matrix(:,:,p));
        percentile=prctile(matrix,threshold,"all");
        for i=1:n_nodes-1
            for j=i:n_nodes
                if matrix(i,j)>percentile
                    significance_mask(i,j)=1;
                    significance_mask(j,i)=1;
                end
            end
        end
    end
end


