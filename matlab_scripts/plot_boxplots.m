function [] = plot_boxplots(metrics,idx_metrics,idx_groups,metrics_labels,patient_labels)
% This function plots the boxplots of the metrics given by the indexes in
% idx_metrics of the groups in idx_groups. For labelling purposes it also
% uses the metrics_labels and the patient_labels.


labels=patient_labels(idx_groups);

n_metrics=length(idx_metrics);
for i=1:n_metrics
    
    % Get data
    metric=idx_metrics(i);
    x=[];
    group=[];
    for g=1:length(idx_groups)
        x=[x metrics{idx_groups(g)}(metric,:)];
        group=[group (g-1)*ones(size(metrics{idx_groups(g)}(metric,:)))];
    end
    
    % Plot
    figure('color','w');
    boxplot(x,group,'Labels',labels)
    title(metrics_labels(metric),'interpreter', 'none','FontSize',15,'FontWeight','normal','FontName','Arial')
    if log10(max(x))-log10(min(x))>2
        set(gca, 'YScale', 'log')
    end
    set(gca,'FontSize',15)
end
end

