function [metrics] = calculate_metrics_v2(mat)
% This function calculates the graph metrics given a connectivity matrix

len_mat=1./mat;                 % conection-length matrix
d_mat= distance_wei(len_mat);   % distance matrix

%calculate metrics
%BC=betweenness_wei(len_mat)';                                   % betweenness centrality
[L,GE]=charpath(d_mat);                                         % characteristic path length and global efficiency
Ci=clustering_coef_wu(weight_conversion(mat, 'normalize'));    % local clustering coefficient
C=mean(Ci);                                                     % global clustering coefficient                            % eigenvector centrality
[Qi, Q] = modularity_und(mat);Qi=Qi';                           % modularity                                   % node strength
T=transitivity_wu(mat);                                         % transitivity
S=smallworldness(mat);                                          %smallworldness

% insert metrics in matrix
metrics=[L GE C Q T S]';
end

