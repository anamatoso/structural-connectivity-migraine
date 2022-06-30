function [] = plot_boxplots_conn(connectomes,idx_metrics,idx_groups,patient_labels,node_labels)

n_nodes=length(node_labels);
n_comparisons=(n_nodes*n_nodes-n_nodes)/2;
labels=patient_labels(idx_groups);

for i=1:length(idx_metrics)
    index1=idx_metrics(i,1);
    index2=idx_metrics(i,2);
    figure('color','w');

    x=[];
    group=[];
    for g=1:length(idx_groups)
        x=[x squeeze(connectomes{g}(index1,index2,:))'];
        group=[group (g-1)*ones(size(squeeze(connectomes{g}(index1,index2,:))'))];
    end
    
    boxplot(x,group,'Labels',labels)
    title(node_labels(index1)+"-"+node_labels(index2),'interpreter', 'none')
%     [h,p]=ttest2(hc_mid,mig_inter,'Alpha',0.05/n_comparisons);
%     if h==1
%         hold on
%         text(2.1,1.01*quantile(mig_inter,0.75), '*','FontSize',15,'Color','black');
%     end
%     text(2.2,1.01*max(max(mig_inter),max(hc_mid)), "p="+num2str(p*n_comparisons),'FontSize',10,'Color','black');
%     hold off
end

end
