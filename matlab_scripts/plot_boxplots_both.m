function [] = plot_boxplots_both(metrics,idx_metrics,ybound,yscale)
idx_groups1=[1 2 5 6];
idx_groups2=[3 4 7 8];

metrics_labels=["Characteristic Path Length" "Global Efficiency" "Clustering Coefficient" "Modularity" "Average Strength" "Transitivity" "Small-worldness"];

patient_labels=["MRtrix-HC-midcycle" "MRtrix-M-interictal" "MRtrix-HC-premenstrual" "MRtrix-M-ictal" "FSL-HC-midcycle" "FSL-M-interictal" "FSL-HC-premenstrual" "FSL-M-ictal"];

nmetrics=length(idx_metrics);
for i=1:nmetrics
    figure('color','w','Position',[200 100 1500 500]); %360 277 619 420
    set(gca,'FontSize',15)
    metric=idx_metrics(i);
    x1=[];
    group1=[];
    for g=1:length(idx_groups1)
        x1=[x1 metrics{idx_groups1(g)}(metric,:)];
        group1=[group1 (g-1)*ones(size(metrics{idx_groups1(g)}(metric,:)))];
    end
    x2=[];
    group2=[];
    for g=1:length(idx_groups2)
        x2=[x2 metrics{idx_groups2(g)}(metric,:)];
        group2=[group2 (g-1)*ones(size(metrics{idx_groups2(g)}(metric,:)))];
    end
    group2=group2+4;
    
    yyaxis left
    boxplot([x1 nan(1,100)],[group1 4*ones(1,25) 5*ones(1,25) 6*ones(1,25) 7*ones(1,25)],"Colors","b")%,'Labels',labels1)
    if yscale
        set(gca, 'YScale', 'log')
    end
    if ~isempty(ybound)
        ylim(ybound)
    end
    ylim([7500 11500])
    yyaxis right
    boxplot([nan(1,100) x2],[0*ones(1,25) 1*ones(1,25) 2*ones(1,25) 3*ones(1,25) group2],"Colors","r",'Labels',patient_labels)
    if yscale
        set(gca, 'YScale', 'log')
    end
    if ~isempty(ybound)
        ylim(ybound)
    end
    title(metrics_labels(metric),'interpreter', 'none','FontSize',20,'FontWeight','normal','FontName','Arial')


    %Image = getframe(gcf);
    %imwrite(Image.cdata, "/Users/ana/Desktop/"+metrics_labels(metric)+".png");
end

end

