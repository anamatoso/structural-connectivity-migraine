function [metrics] = calculate_metrics(mat,version)
% This function calculates the graph metrics given a connectivity matrix

if version==1
    
    len_mat=1./mat;                 % conection-length matrix
    d_mat= distance_wei(len_mat);   % distance matrix
    
    %calculate metrics
    BC=betweenness_wei(len_mat)';                               % betweenness centrality 1-116
    [L,GE]=charpath(d_mat,0,0);                                     % characteristic path length and global efficiency 117,118
    Ci=clustering_coef_wu(weight_conversion(mat, 'normalize'))';% local clustering coefficient 120-235
    C=mean(Ci);                                                 % global clustering coefficient 119
    EC=eigenvector_centrality_und(mat)';                        % eigenvector centrality 236-351
    [Qi, Q]=modularity_und(mat);Qi=Qi';                         % modularity 352-467, 468
    RC=rich_club_wu_norm(mat,115);                              % rich club coefficient 469-583
    strength=strengths_und(mat);                                % node strength 584-699
    mean_strength=mean(strength);                               % mean strength 700
    T=transitivity_wu(mat);                                     % transitivity 701
    S=smallworldness(mat);                                      % smallworldness 702
    A=assortativity_wei(mat,0);                                 % assortivity 703
    
    % insert metrics in matrix
    metrics=[BC L GE C Ci EC Qi Q RC strength mean_strength T S A]';
    
    
elseif version==2
    len_mat=1./mat;                 % conection-length matrix
    d_mat= distance_wei(len_mat);   % distance matrix
    
    %calculate metrics
    [L,GE]=charpath(d_mat,0,0);                                     % characteristic path length and global efficiency 1,2
    Ci=clustering_coef_wu(weight_conversion(mat, 'normalize'))';    % local clustering coefficient 
    C=mean(Ci);                                                     % global clustering coefficient 3    
    [~, Q]=modularity_und(mat);                                          % modularity 4
    strength=strengths_und(mat);                                    % node strength 5-120
    mean_strength=mean(strength);                                   % mean strength 121
    T=transitivity_wu(mat);                                         % transitivity 122
    S=smallworldness(mat);                                          % smallworldness 123
    A=assortativity_wei(mat,0);                                     % assortivity 124
    
    % insert metrics in matrix
    metrics=[L GE C Q strength mean_strength T S A]';


elseif version==3
    len_mat=1./mat;                 % conection-length matrix
    d_mat= distance_wei(len_mat);   % distance matrix
    
    %calculate metrics
    [L,GE]=charpath(d_mat,0,0);                                     % characteristic path length and global efficiency 1,2
    Ci=clustering_coef_wu(weight_conversion(mat, 'normalize'));    % local clustering coefficient 
    C=mean(Ci);                                                     % global clustering coefficient 3    
    [~, Q]=modularity_und(mat);                                          % modularity 4
    strength=strengths_und(mat);                                    % node strength 5-120
    mean_strength=mean(strength);                                   % mean strength 121
    T=transitivity_wu(mat);                                         % transitivity 122
    S=smallworldness(mat);                                          % smallworldness 123
    A=assortativity_wei(mat,0);                                     % assortivity 124
    
    % insert metrics in matrix
    metrics=[L GE C Q mean_strength T S A]';
end
end

