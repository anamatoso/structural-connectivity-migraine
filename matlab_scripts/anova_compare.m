function [ANOVA_results] = anova_compare(metrics,metrics_labels,version_metrics,n_nodes,varargin)
% Given the metrics and the respective labels, this function calculates
% which metrics are significantly different. If varargin is set to "True"
% then it shows all metrics (instead of only the significant ones).
n_metrics=length(metrics_labels);
compare_anova=zeros(1,5);

if version_metrics==3
    correction=1;
else
    correction=n_nodes;%n_metrics;
end

if ~isempty(varargin)
    if varargin{1}=="True"
        threshold=inf;
    else
        threshold=0.05/correction;
    end
end

for m=1:n_metrics
    hc_mid=metrics{1}(m,:);
    mig_inter=metrics{2}(m,:);
    hc_pre=metrics{3}(m,:);
    mig_ict=metrics{4}(m,:);
    
    x = [hc_mid mig_inter];% hc_pre mig_ict];
    g = [zeros(1,length(hc_mid)),ones(1,length(mig_inter))];%,2.*ones(1,length(hc_pre)),3.*ones(1,length(mig_ict))];
    
    if any(isnan(x)) || all(x==0)
        compare_anova=[compare_anova;m 1 2 correction 0];
        continue
    end
    
    [~, ~, stats] = anova1(x,g,'off');
    c=multcompare(stats,'display','off','alpha',0.05/correction,'ctype','bonferroni');
    
    for idx_p=1:length(c(:,end))
        if c(idx_p,end)<threshold
            switch c(idx_p,4)>0 % difference is positive?
                case 1
                    compare_anova=[compare_anova;m c(idx_p,1) c(idx_p,2) c(idx_p,end) 1];
                case 0
                    compare_anova=[compare_anova;m c(idx_p,1) c(idx_p,2) c(idx_p,end) -1];
            end
        end
    end
    
end

compare_anova=compare_anova(2:end,:);
table=array2table(compare_anova, "VariableNames", ["Metric index","Group 1", "Group 2", "P-value", "Difference"]);
metrics_names=metrics_labels(compare_anova(:,1))';
t_names=array2table(metrics_names, "VariableNames", ["Metric Name"]);
ANOVA_results = [t_names table];
end

