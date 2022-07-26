function [ANOVA_results] = anova_compare(metrics,metrics_labels,version_metrics,n_nodes)
% Given the metrics and the respective labels, this function calculates
% which metrics are significantly different. If varargin is set to "True"
% then it shows all metrics (instead of only the significant ones).
n_metrics=length(metrics_labels);
compare_anova=zeros(1,5);

if version_metrics==2
    correction=1;
else
    correction=n_nodes;%n_metrics;
end

for m=1:n_metrics
    data1=metrics{1}(m,:);
    data2=metrics{2}(m,:);
    data3=metrics{3}(m,:);
    data4=metrics{4}(m,:);

    x = [data1 data2 data3 data4];
    g = [ones(1,length(data1)),2*ones(1,length(data2)),...
        3*ones(1,length(data2)),4*ones(1,length(data2))];

    if any(isnan(x))
        compare_anova=[compare_anova;m 1 2 correction 0];
        continue
    elseif all(x==0)
        compare_ttest=[compare_ttest;m g1 g2 0 0 0];
        continue
    end
    
    allnorm=(1-kstest(data1))*(1-kstest(data2))*(1-kstest(data3))*(1-kstest(data4));
    
    if kstest(data1) && 
    [~, ~, stats] = anova1(x,g,'off');
    
    
    c=multcompare(stats,'display','off','alpha',0.05/correction,'ctype','bonferroni');

    for idx_p=1:length(c(:,end))
        switch c(idx_p,4)>0 % difference is positive?
            case 1
                compare_anova=[compare_anova;m c(idx_p,1) c(idx_p,2) c(idx_p,end) 1];
            case 0
                compare_anova=[compare_anova;m c(idx_p,1) c(idx_p,2) c(idx_p,end) -1];
        end
    end
end

compare_anova=compare_anova(2:end,:);
table=array2table(compare_anova, "VariableNames", ["Metric index","Group 1", "Group 2", "P-value", "Difference"]);
metrics_names=metrics_labels(compare_anova(:,1))';
t_names=array2table(metrics_names, "VariableNames", ["Metric Name"]);
ANOVA_results = [t_names table];
end

