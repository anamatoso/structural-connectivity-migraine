function [ttest_results] = ttest_compare(metrics,metrics_labels,version_metrics,n_nodes,comparisons)
% Given the metrics and the respective labels, this function calculates
% which metrics are significantly different.

n_metrics=length(metrics_labels);
compare_ttest=zeros(1,6);

if version_metrics==2
    correction=1;
else
    correction=n_nodes;
end
[ncomp,~]=size(comparisons);
for c=1:ncomp
    g1=comparisons(c,1);
    g2=comparisons(c,2);
    for m=1:n_metrics
        data1=metrics{g1}(m,:);
        data2=metrics{g2}(m,:);
        
        x = [data1 data2];
        g = [zeros(1,length(data1)),ones(1,length(data2))];
        
        if any(isnan(x)) || all(x==0)
            compare_ttest=[compare_ttest;m g1 g2 1 correction 0];
            continue
        end
        
        [~,p] = ttest2(data1,data2);
        switch mean(data2)>mean(data1)
            case 1
                compare_ttest=[compare_ttest;m g1 g2 p p*correction 1];
            case 0
                compare_ttest=[compare_ttest;m g1 g2 p p*correction -1];
        end
    end
end

compare_ttest=compare_ttest(2:end,:);
table=array2table(compare_ttest, "VariableNames", ["Metric index","Group 1", "Group 2", "P-value","P-value (corrected)", "Difference"]);
metrics_names=metrics_labels(compare_ttest(:,1))';
t_names=array2table(metrics_names, "VariableNames", ["Metric Name"]);
ttest_results = [t_names table];

%[fdr,q] = mafdr(table2array(ttest_results(:,5)));
%fdr_table=array2table([fdr q], "VariableNames", ["FDR", "Q"]);

ttest_results = [t_names table];
if version_metrics==3
    % Calculate FDR for nodal metrics
    pvals=reshape(table2array(ttest_results(:,5)),116,[]);
    [~,nmetrics]=size(pvals);
    fdr=zeros(size(pvals));
    for i=1:nmetrics
        [~, q] = mafdr(pvals(:,i));
        fdr(:,i)=q;
    end
    
    fdr=reshape(fdr,[],1);
    fdr_table=array2table(fdr, "VariableNames", ["q-value"]);
    ttest_results = [t_names table fdr_table];
end

end

