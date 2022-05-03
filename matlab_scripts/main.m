%% Load data from matrices
% directory where connectivity matrices are
dir='/Users/ana/Documents/Ana/universidade/Tese/Code/matlab_scripts/mat_data';

% Controls midcyle
HC_midcycle_mrtrix=load_data(dir,'*midcycle*mrtrix*'); %116 x 116 x n_people
HC_midcycle_fsl=load_data(dir,'*midcycle*fsl*'); 

% Patients interictal
M_interictal_mrtrix=load_data(dir,'*interictal*mrtrix*');
M_interictal_fsl=load_data(dir,'*interictal*fsl*');

connectomes={HC_midcycle_mrtrix HC_midcycle_fsl M_interictal_mrtrix M_interictal_fsl};

% Calculate people per situation
n_people=zeros(1,length(connectomes));
for i=1:length(connectomes)
    conmat=connectomes{i};
    s=size(conmat);
    n_people(i)=s(end);
end

clear dir s conmat i
%% Test for spurious connections

significance_mask=zeros(116,116,length(connectomes));
for i=1:length(connectomes)
    significance_mask(:,:,i) = signtest_mask(connectomes(i));    
end
% imagesc(significance_mask(i)); colormap jet
%% Calculate metrics

for i=1:length(connectomes)
    conmats=connectomes{i};
    for p=1:n_people(i)
        mat=conmats(p);                 % connectivity matrix
        len_mat=1./mat;                 % conection-length matrix
        d_mat= distance_wei(len_mat);   % distance matrix
        
        %calculate metrics
        BC=betweenness_wei(len_mat);                                % betweenness centrality
        [L,GE]=charpath(d_mat);                                         % characteristic path length and global efficiency
        Ci=clustering_coef_wu(weight_conversion(mat, 'normalize')); % local clustering coefficient
        C=mean(Ci);                                                 % global clustering coefficient
        EC=eigenvector_centrality_und(mat);                         % eigenvector centrality
        [Qi, Q] = modularity_und(mat);                              % modularity
        RC=rich_club_wu(mat);                                          % rich club coefficient
        %degree=strengths_und(mat);                                  % node strength
        T=transitivity_wu(mat);                                     % transitivity
        S=smallworldness(C,L,mat);                                  %smallworldness
    
        % insert metrics in matrix
        metrics{i,p}=[BC L GE C EC Q RC T S];     
    end
end

%% Analysis of results


