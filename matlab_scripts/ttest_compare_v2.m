function [ttest_results] = ttest_compare_v2(metrics,metrics_labels,version_metrics,n_nodes,comparisons)
% Given the metrics and the respective labels, this function calculates
% which metrics are significantly different.
cycle=[1 5;3 7;2 6;4 8];

if all(comparisons==cycle)
    paired=true;
else 
    paired=false;
end

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
        %g = [zeros(1,length(data1)),ones(1,length(data2))];
        
        if any(isnan(x)) || all(x==0)
            compare_ttest=[compare_ttest;m g1 g2 1 correction 0];
            continue
        end
        
        if ~kstest(data1) && ~kstest(data2) && ~vartest2(data1,data2)   % normal with equal variances
            if paired
                [~,p] = ttest(data1,data2);
            else
                [~,p] = ttest2(data1,data2);
            end
            dif=mean(data2)-mean(data1);
        
        elseif ~kstest(data1) && ~kstest(data2) && vartest2(data1,data2)% normal with unequal variances
            [~,p] = ttest2(data1,data2,'vartype', 'unequal');
            dif=mean(data2)-mean(data1);
        else                                                            % not normal with unequal variances
            %if paired
             %   p = signrank(data1,data2);
            %else
                p = ranksum(data1,data2);
            %end
            dif=median(data2)-median(data1);
        end
        
        switch dif>0
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

if version_metrics==3
    % Calculate FDR for nodal metrics
    pvals=reshape(table2array(ttest_results(:,5)),n_nodes,[]);
    [~,nmetrics]=size(pvals);
    fdr=zeros(size(pvals));
    for i=1:nmetrics
        [q] = mafdr(pvals(:,i),'BHFDR',true);
        %[fdr2_v,q2_v] = mafdr(pvals(:,i));
        fdr(:,i)=q;
        %fdr2(:,i)=fdr2_v;
        %q2(:,i)=q2_v;
    end
    
    fdr=reshape(fdr,[],1);
    %fdr2=reshape(fdr2,[],1);
    %q2=reshape(q2,[],1);
    %fdr_table=array2table([fdr fdr2 q2], "VariableNames", ["q-value(BH)" "FDR" "qvalue(S)"]);
    fdr_table=array2table(fdr, "VariableNames", ["q-value"]);
    ttest_results = [t_names table fdr_table];
end

end

