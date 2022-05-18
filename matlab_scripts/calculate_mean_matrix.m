function [matrices] = calculate_mean_matrix(connectomes)
% this function calculates the mean connectivity matrix of each group by
% averaging among the subjects of each group.

n_conditions=length(connectomes);
matrices=zeros(116,116,n_conditions);
for condition=1:n_conditions
   m=connectomes{condition};
   matrices(:,:,condition)=mean(m,3);
end
end

