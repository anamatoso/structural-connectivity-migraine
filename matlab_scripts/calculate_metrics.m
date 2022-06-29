function [metrics] = calculate_metrics(mat,version)
% This function calculates the graph metrics given a connectivity matrix and a given choice of metrics

n_nodes=length(mat);
len_mat=1./mat; % conection-length matrix
if version==1
  
    % metrics
    BC=betweenness_wei(len_mat)'/((n_nodes-1)*(n_nodes-2));     % betweenness centrality 1-116
    Ci=clustering_coef_wu(weight_conversion(mat, 'normalize'))';% local clustering coefficient 117-232
    EC=eigenvector_centrality_und(mat)';                        % eigenvector centrality 233-348
    RC=rich_club_wu_norm(mat,n_nodes-1);                        % rich club coefficient 249-463
    strength=strengths_und(mat);                                % node strength 464-579
    
    % insert metrics in matrix
    metrics=[BC Ci EC RC strength]';
    
elseif version==2
    d_mat= distance_wei(len_mat);   % distance matrix
    
    %calculate metrics
    [L,GE]=charpath(d_mat,0,0);                                     % characteristic path length and global efficiency 1,2
    Ci=clustering_coef_wu(weight_conversion(mat, 'normalize'));     % local clustering coefficient 
    C=mean(Ci);                                                     % global clustering coefficient 3    
    [~, Q]=modularity_und(mat);                                     % modularity 4
    strength=strengths_und(mat);                                    % node strength 
    mean_strength=mean(strength);                                   % mean strength 5
    T=transitivity_wu(weight_conversion(mat, 'normalize'));         % transitivity 6
    S=smallworldness(mat);                                          % smallworldness 7
    %A=assortativity_wei(mat,0);                                     % assortivity 8
    
    % insert metrics in matrix
    metrics=[L GE C Q mean_strength T S]';


end

