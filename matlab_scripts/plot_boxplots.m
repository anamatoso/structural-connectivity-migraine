function [] = plot_boxplots(metrics,idx,metrics_labels,version_metrics,n_nodes)
n_metrics=length(metrics_labels);

if version_metrics==3
    correction=1;
else
    correction=n_nodes;
end


for i=1:length(idx)
    index=idx(i);
    figure('color','w');
    
    hc_mid=metrics{1}(index,:);
    mig_inter=metrics{2}(index,:);
    
    x = [hc_mid mig_inter];
    g = [zeros(1,length(hc_mid)),ones(1,length(mig_inter))];
    
    boxplot(x,g,'Labels',{'HC','M'})
    title(metrics_labels(index),'interpreter', 'none','FontSize',15,'FontWeight','normal','FontName','Arial')
    [h,p]=ttest2(hc_mid,mig_inter,'Alpha',0.05/correction);
    if h==1
        hold on
        %text(2.1,1.01*quantile(mig_inter,0.75), '*','FontSize',15,'Color','black');
    end
    text(2.2,1.01*max(max(mig_inter),max(hc_mid)), "p="+num2str(p*correction,2),'FontSize',15,'Color','black');
    hold off
end

end

