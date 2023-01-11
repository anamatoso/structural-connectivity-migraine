function [metrics] = calculate_metrics_v2(mat,v)
% This function calculates the graph metrics given a connectivity matrix and a given choice of metrics
n_nodes=length(mat);
len_mat=1./mat; % conection-length matrix
d_mat= distance_wei(len_mat);   % distance matrix

randmat=random_matrix(mat);
len_randmat=1./randmat; % conection-length matrix
d_randmat= distance_wei(len_randmat);   % distance matrix

if v==2
    %calculate metrics
    [L,GE]=charpath(d_mat,0,0);                                     % characteristic path length and global efficiency 1,2
    [Lrand,GErand]=charpath(d_randmat,0,0);  
    
    Ci=clustering_coef_wu(weight_conversion(mat, 'normalize'));     % local clustering coefficient 
    C=mean(Ci);                                                     % global clustering coefficient 3    
    
    [~, Q]=modularity_und(mat);                                     % modularity 4
    
    strength=strengths_und(mat);                                    % node strength 
    mean_strength=mean(strength);                                   % mean strength 5
    strengthrand=strengths_und(randmat);
    mean_strengthrand=mean(strengthrand);
    
    T=transitivity_wu(weight_conversion(mat, 'normalize'));         % transitivity 6
    S=smallworldness3(mat);                                         % smallworldness 7
    
    
    % insert metrics in matrix
    metrics=[L/Lrand GE/GErand C Q mean_strength/mean_strengthrand T S]';
else
    % metrics
    BC=betweenness_wei(len_mat)'/((n_nodes-1)*(n_nodes-2));         % betweenness centrality 1-116
    Ci=clustering_coef_wu(weight_conversion(mat, 'normalize'))';    % local clustering coefficient 117-232
    EC=eigenvector_centrality_und(mat)';                            % eigenvector centrality 233-348
    strength=strengths_und(mat);                                    % node strength 464-579
    strengthrand=strengths_und(randmat);

    % insert metrics in matrix
    metrics=[BC Ci EC strength./strengthrand]';

end



end

