function [ttest_results] = ttest_compare_v2(metrics,metrics_labels,n_nodes,comparisons)
% Given the metrics and the respective labels, this function calculates the
% pvalues of the comparison of all metrics between each pair of groups in
% comparisons. It outputs a table with the following headers: Metric Name,
% Group 1, Group 2, P-value, P-value (corrected) and Difference. If the
% metrics are nodal then it also has another column with the FDR corrected
% values.

n_metrics=numel(metrics_labels);
[n_comp,~]=size(comparisons);
compare_ttest=zeros(n_metrics*n_comp,6);

% check whether a correction must be done
if n_metrics==7 %(aka version metrics=2)
    correction=1;
else
    correction=n_nodes;
end


idx=1;% index of table
for c=1:n_comp

    % Select groups to compare
    g1=comparisons(c,1);
    g2=comparisons(c,2);

    % Go through all the metrics
    for m=1:n_metrics

        % Select the metrics of those groups
        data1=metrics{g1}(m,:);
        data2=metrics{g2}(m,:);
        x = [data1 data2];
        
        % If there is there is NaN or if it is all==0, then it is not
        % interesting, hence p is maximum and difference is set to 0
        if any(isnan(x)) || all(x==0)
            compare_ttest(idx,:)=[m g1 g2 1 correction 0];
            idx=idx+1;
            continue
        end
        
        % Select the most appropriate statistical test
        
        % Option 1 - data is normal and has equal variances
        if ~kstest(data1) && ~kstest(data2) && ~vartest2(data1,data2)   
            [~,p] = ttest2(data1,data2);
            dif=mean(data2)-mean(data1);
        
        % Option 2 - data is normal and has unequal variances
        elseif ~kstest(data1) && ~kstest(data2) && vartest2(data1,data2)
            [~,p] = ttest2(data1,data2,'vartype', 'unequal');
            dif=mean(data2)-mean(data1);
        
        % Option 3 - data is not normal and has unequal variances
        else                                                            
            p = ranksum(data1,data2);
            dif=median(data2)-median(data1);
        end
        
        % Determine if the trend is increased or decreased in the second
        % group
        switch dif>0
            case 1
                compare_ttest(idx,:)=[m g1 g2 p p*correction 1];
            case 0
                compare_ttest(idx,:)=[m g1 g2 p p*correction -1];
        end
        idx=idx+1;
    end
end

% Create tables
table1=array2table(compare_ttest(:,2:end), "VariableNames", ["Group 1", "Group 2", "P-value","P-value (corrected)", "Difference"]);
metrics_names=array2table(metrics_labels(compare_ttest(:,1))', "VariableNames", "Metric Name");
ttest_results = [metrics_names table1]; % Join tables


% If we are dealing with nodal metrics, perform FDR correction and add
% column with q-values
if n_metrics~=7
    
    % Calculate FDR for nodal metrics
    pvals=reshape(table2array(ttest_results(:,5)),n_nodes,[]);
    [~,nmetrics]=size(pvals);
    FDR=zeros(size(pvals));
    for i=1:nmetrics
        q = mafdr(pvals(:,i),'BHFDR',true);
        FDR(:,i)=q;
    end
    
    % Insert column in table
    FDR=reshape(FDR,[],1);
    FDR_table=array2table(FDR, "VariableNames", "q-value");
    ttest_results = [ttest_results FDR_table];
end
end

