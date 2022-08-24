function [] = plot_boxplots(metrics,idx_metrics,idx_groups,metrics_labels,patient_labels,version_metrics,n_nodes)

if version_metrics==2
    correction=1;
else
    correction=n_nodes;
end

labels=patient_labels(idx_groups);

nmetrics=length(idx_metrics);
for i=1:nmetrics
    figure('color','w');
    metric=idx_metrics(i);
    x=[];
    group=[];
    for g=1:length(idx_groups)
        x=[x metrics{idx_groups(g)}(metric,:)];
        group=[group (g-1)*ones(size(metrics{idx_groups(g)}(metric,:)))];
    end
    
    

    boxplot(x,group,'Labels',labels)
    title(metrics_labels(metric),'interpreter', 'none','FontSize',15,'FontWeight','normal','FontName','Arial')
%     for g1=1:length(idx_groups)-1
%         for g2=g1+1:length(idx_groups)
%             [h,~]=ttest2(metrics{idx_groups(g1)}(metric,:),metrics{idx_groups(g2)}(metric,:),'Alpha',0.05/correction);
%             if h==1
%                 yt = max(max(metrics{idx_groups(g1)}(metric,:)),max(metrics{idx_groups(g2)}(metric,:)));%get(gca, 'YTick');
%                 axis([xlim    0  max(yt)*1.2])
%                 xt = get(gca, 'XTick');
%                 hold on
%                 plot(xt([g1+0.1 g2-0.1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([g1 g2])), max(yt)*1.15, '*k')
%                 hold off
%                 
%             end
%         end
%     end
    if log(max(x))-log(min(x))>2
        set(gca, 'YScale', 'log')
    end
    
    
    
    
    %
    %      if h==1
    %          hold on
    %          text(2.1,1.01*quantile(g2,0.75), '*','FontSize',15,'Color','black');
    %      end
    %text(2.2,1.01*max(max(mig_inter),max(hc_mid)), "p="+num2str(p*correction,2),'FontSize',15,'Color','black');
    %     hold off
    %ax = gca;
    %ax.YAxis.Scale ="log";
end

end

