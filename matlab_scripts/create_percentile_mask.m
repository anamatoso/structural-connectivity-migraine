function [significance_mask] = create_percentile_mask(group_matrix,threshold)
% This function creates a percentile mask where 1 means that the connection is above the threshold and 0 otherwise.

[n_nodes,~,n_people]=size(group_matrix);
significance_mask=zeros(size(group_matrix));

for p=1:n_people
    matrix=squeeze(group_matrix(:,:,p));
    percentile=prctile(matrix,threshold,"all");
    for i=1:n_nodes-1
        for j=i:n_nodes
            if matrix(i,j)>percentile
                significance_mask(i,j,p)=1;
                significance_mask(j,i,p)=1;
            end
        end
    end
end



