function [] = plot_boxplots(metrics,idx,metrics_labels)
n_metrics=length(metrics_labels);
for i=1:length(idx)
    index=idx(i);
    figure;
    
    hc_mid=metrics{1}(index,:);
    mig_inter=metrics{2}(index,:);
    
    x = [hc_mid mig_inter];
    g = [zeros(1,length(hc_mid)),ones(1,length(mig_inter))];
    
    boxplot(x,g,'Labels',{'HC','M'})
    title(metrics_labels(index),'interpreter', 'none')
    [h,p]=ttest2(hc_mid,mig_inter,'Alpha',0.05/n_metrics);
    if h==1
        hold on
        text(2.1,1.01*quantile(mig_inter,0.75), '*','FontSize',15,'Color','black');
    end
    text(2.2,max(max(mig_inter),max(hc_mid)), "p="+num2str(p*n_metrics),'FontSize',10,'Color','black');

    hold off
end

end

