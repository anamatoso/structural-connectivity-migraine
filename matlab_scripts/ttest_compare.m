function [ANOVA_results] = ttest_compare(metrics,metrics_labels,version_metrics,n_nodes,comparisons,pvalue)
% Given the metrics and the respective labels, this function calculates
% which metrics are significantly different. 

n_metrics=length(metrics_labels);
compare_ttest=zeros(1,5);

if version_metrics==2
    correction=1;
else
    correction=n_nodes;%n_metrics;
end
[ncomp,~]=size(comparisons);
for c=1:ncomp
    for m=1:n_metrics
        g1=comparisons(c,1);
        g2=comparisons(c,2);
        
        data1=metrics{g1}(m,:);
        data2=metrics{g2}(m,:);  
        
        x = [data1 data2];% hc_pre mig_ict];
        g = [zeros(1,length(data1)),ones(1,length(data2))];
        
        if any(isnan(x)) || all(x==0)
            compare_ttest=[compare_ttest;m g1 g2 correction 0];
            continue
        end
        
        [~,p] = ttest2(data1,data2);
        if pvalue=="pcorrected"
            p=p*correction;
        end
        switch mean(data2)>mean(data1)
            case 1
                compare_ttest=[compare_ttest;m g1 g2 p 1];
            case 0
                compare_ttest=[compare_ttest;m g1 g2 p -1];
        end  
    end
end
if pvalue=="pcorrected"
    pvalue="P-value (corrected)";
else
    pvalue="P-value";
end
compare_ttest=compare_ttest(2:end,:);
table=array2table(compare_ttest, "VariableNames", ["Metric index","Group 1", "Group 2", pvalue, "Difference"]);
metrics_names=metrics_labels(compare_ttest(:,1))';
t_names=array2table(metrics_names, "VariableNames", ["Metric Name"]);
ANOVA_results = [t_names table];
end

