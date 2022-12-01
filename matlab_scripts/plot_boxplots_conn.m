function [] = plot_boxplots_conn(connectomes,connectivity_pairs,idx_groups,patient_labels,node_labels)
% This function plots the boxplots of the connectivity given by the pairs
% of indexes in connectivity_pairs of the groups in idx_groups. For
% labelling purposes it also uses the node_labels and the
% patient_labels.


labels=patient_labels(idx_groups);

for i=1:length(connectivity_pairs)
    
    % Get data
    index1=connectivity_pairs(i,1);
    index2=connectivity_pairs(i,2);
    
    x=[];
    group=[];
    for g=1:length(idx_groups)
        x=[x squeeze(connectomes{g}(index1,index2,:))'];
        group=[group (g-1)*ones(size(squeeze(connectomes{g}(index1,index2,:))'))];
    end
    
    % Plot
    figure('color','w');
    boxplot(x,group,'Labels',labels)
    title(node_labels(index1)+"-"+node_labels(index2),'interpreter', 'none')
    set(gca,'FontSize',15)
end
end
