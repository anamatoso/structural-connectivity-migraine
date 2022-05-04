function [metrics] = calculate_metrics(mat)
% This function calculates the graph metrics given a connectivity matrix

len_mat=1./mat;                 % conection-length matrix
d_mat= distance_wei(len_mat);   % distance matrix

%calculate metrics
BC=betweenness_wei(len_mat)';                                   % betweenness centrality
[L,GE]=charpath(d_mat);                                         % characteristic path length and global efficiency
Ci=clustering_coef_wu(weight_conversion(mat, 'normalize'))';    % local clustering coefficient
C=mean(Ci);                                                     % global clustering coefficient
EC=eigenvector_centrality_und(mat)';                            % eigenvector centrality
[Qi, Q] = modularity_und(mat);Qi=Qi';                           % modularity
RC=rich_club_wu(mat);                                           % rich club coefficient
%degree=strengths_und(mat);                                     % node strength
T=transitivity_wu(mat);                                         % transitivity
S=smallworldness(mat);                                          %smallworldness

% insert metrics in matrix
metrics=[BC L GE C Ci EC Q Qi RC T S]';
end

