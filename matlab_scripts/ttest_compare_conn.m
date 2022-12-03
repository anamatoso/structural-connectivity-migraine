function [ttest_results_conn] = ttest_compare_conn(connectomes,node_labels,comparisons)
% This function performs the statistical tests and saves the properties in
% a table with the following headers: Node 1, Node 2, Group 1, Group 2,
% P-value,P-value (corrected), Difference and q- values.

n_nodes=length(node_labels);
n_pairs=(n_nodes*n_nodes-n_nodes)/2;
[n_comp,~]=size(comparisons);
compare_ttest=zeros(n_pairs*n_comp,7);

idx=1;% index of table
for c=1:n_comp

    % Select groups to compare
    g1=comparisons(c,1);
    g2=comparisons(c,2);
    
    % Iterate through all pairs
    for m=1:n_nodes-1
        for n=m+1:n_nodes
            
            % Select the connnectivity of the pairs of those groups
            data1=squeeze(connectomes{g1}(m,n,:))';
            data2=squeeze(connectomes{g2}(m,n,:))';
            
            % If there is there is NaN or if it is all==0, then it is not
            % interesting, hence p is maximum and difference is set to 0
            x = [data1 data2];
            if any(isnan(x)) || all(x==0)
                compare_ttest(idx,:)=[m n g1 g2 1 n_pairs 0];
                idx=idx+1;
                continue
            end
            
            % Perform ttest and see if it is increased or decreased in
            % the second group
            [~,p] = ttest2(data1,data2);
            switch mean(data2)>mean(data1)
                case 1
                    compare_ttest(idx,:)=[m n g1 g2 p p*n_pairs 1];
                case 0
                    compare_ttest(idx,:)=[m n g1 g2 p p*n_pairs -1];
            end
            idx=idx+1;
        end
    end
end

% Create tables
table=array2table(compare_ttest, "VariableNames", ["Node 1","Node 2","Group 1", "Group 2", "P-value","P-value (corrected)","Difference"]);
node_names=[node_labels(compare_ttest(:,1));node_labels(compare_ttest(:,2))]';
node_names_table=array2table(node_names, "VariableNames", ["Node 1 Name","Node 2 Name"]);
ttest_results_conn = [node_names_table table];


% Calculate FDR for each matrix
pvals=reshape(table2array(ttest_results_conn(:,7)),[],n_comp);
[~,ncomparisons]=size(pvals);
FDR=zeros(size(pvals));
for i=1:ncomparisons
    q = mafdr(pvals(:,i),'BHFDR',true);
    FDR(:,i)=q;
end

% Insert column in table
FDR=reshape(FDR,[],1);
FDR_table=array2table(FDR, "VariableNames", "q-value");
ttest_results_conn = [ttest_results_conn FDR_table];

end

