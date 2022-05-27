function [] = plot_boxplots_conn(connectomes,idx,node_labels)

n_nodes=length(node_labels);
n_comparisons=(n_nodes*n_nodes-n_nodes)/2;

for i=1:length(idx)
    index1=idx(i,1);
    index2=idx(i,2);
    figure;
    
    hc_mid=squeeze(connectomes{1}(index1,index2,:))';
    mig_inter=squeeze(connectomes{2}(index1,index2,:))';
    
    x = [hc_mid mig_inter];
    g = [zeros(1,length(hc_mid)),ones(1,length(mig_inter))];
    
    boxplot(x,g,'Labels',{'HC','M'})
    title(node_labels(index1)+"-"+node_labels(index2),'interpreter', 'none')
    
    if ttest2(hc_mid,mig_inter,'Alpha',0.05/n_comparisons)==1
        hold on
        text(2.1,1.01*quantile(mig_inter,0.75), '*','FontSize',15,'Color','black');
    end
    hold off
end

end
