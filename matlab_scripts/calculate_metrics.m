function [metrics] = calculate_metrics(mat)
% This function calculates the graph metrics given a connectivity matrix

len_mat=1./mat;                 % conection-length matrix
d_mat= distance_wei(len_mat);   % distance matrix

%calculate metrics
BC=betweenness_wei(len_mat)';                                   % betweenness centrality 1-116
[L,GE]=charpath(d_mat);                                         % characteristic path length and global efficiency 117,118
Ci=clustering_coef_wu(weight_conversion(mat, 'normalize'))';    % local clustering coefficient 120-235
C=mean(Ci);                                                     % global clustering coefficient 119
EC=eigenvector_centrality_und(mat)';                            % eigenvector centrality 236-351
[Qi, Q]=modularity_und(mat);Qi=Qi';                             % modularity 352-467, 468
RC=rich_club_wu_norm(mat,115);                                  % rich club coefficient 469-583
degree=strengths_und(mat);                                      % node strength 584-699
T=transitivity_wu(mat);                                         % transitivity 700
S=smallworldness(mat);                                          % smallworldness 701
A=assortativity_wei(mat,0);
% insert metrics in matrix
metrics=[BC L GE C Ci EC Qi Q RC degree T S A]';
end

