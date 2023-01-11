function [metrics] = calculate_metrics(mat,version)
%CALCULATE_METRICS   Calculates the connectivity metrics of a given matrix.
%   metrics = CALCULATE_METRICS(mat,version) is the vector with the connectivity metrics of mat. 
%   version states which metrics you want to calculate:
%       If version==1 then it calculates the rich club coefficient;
%       If version==2 then it calculates the global metrics;
%       If version==3 then it calculates the nodal metrics.

n_nodes=length(mat);
len_mat=1./mat; % conection-length matrix

if version==1
  
    % metric
    RC=rich_club_wu_norm(mat,n_nodes-1);                            % rich club coefficient
    
    % insert metrics in matrix
    metrics=RC';
    
elseif version==2
    d_mat= distance_wei(len_mat);   % distance matrix
    
    %calculate metrics
    [L,GE]=charpath(d_mat,0,0);                                     % characteristic path length and global efficiency
    Ci=clustering_coef_wu(weight_conversion(mat, 'normalize'));     % local clustering coefficient 
    C=mean(Ci);                                                     % global clustering coefficient 
    [~, Q]=modularity_und(mat);                                     % modularity
    strength=strengths_und(mat);                                    % node strength 
    mean_strength=mean(strength);                                   % mean strength
    T=transitivity_wu(weight_conversion(mat, 'normalize'));         % transitivity
    %S=smallworldness2(mat,20000);                                  % smallworldness
    S=smallworldness3(mat);                                         % smallworldness
    %A=assortativity_wei(mat,0);                                    % assortivity
    
    % insert metrics in matrix
    metrics=[L GE C Q mean_strength T S]';

elseif version==3
  
    % metrics
    BC=betweenness_wei(len_mat)'/((n_nodes-1)*(n_nodes-2));         % betweenness centrality
    Ci=clustering_coef_wu(weight_conversion(mat, 'normalize'))';    % local clustering coefficient
    EC=eigenvector_centrality_und(mat)';                            % eigenvector centrality
    strength=strengths_und(mat);                                    % node strength
    
    % insert metrics in matrix
    metrics=[BC Ci EC strength]';

end

