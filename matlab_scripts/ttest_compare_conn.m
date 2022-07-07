function [ttest_results_conn] = ttest_compare_conn(connectomes,node_labels,comparisons)

compare_ttest=zeros(1,7);
n_nodes=length(node_labels);
n_comparisons=(n_nodes*n_nodes-n_nodes)/2;

[ncomp,~]=size(comparisons);
for c=1:ncomp
    g1=comparisons(c,1);
    g2=comparisons(c,2);
    
    for m=1:n_nodes-1
        for n=m+1:n_nodes
            
            data1=squeeze(connectomes{g1}(m,n,:))';
            data2=squeeze(connectomes{g2}(m,n,:))';
            
            x = [data1 data2];
            g = [zeros(1,length(data1)),ones(1,length(data2))];
            
            if any(isnan(x)) || all(x==0)
                compare_ttest=[compare_ttest;m n g1 g2 1 n_comparisons 0];
                continue
            end
            
            
            [~,p] = ttest2(data1,data2);
            switch mean(data2)>mean(data1)
                case 1
                    compare_ttest=[compare_ttest;m n g1 g2 p p*n_comparisons 1];
                case 0
                    compare_ttest=[compare_ttest;m n g1 g2 p p*n_comparisons -1];
            end
            
        end
    end
end
compare_ttest=compare_ttest(2:end,:);
table=array2table(compare_ttest, "VariableNames", ["Node 1","Node 2","Group 1", "Group 2", "P-value","P-value (corrected)","Difference"]);
node_names=[node_labels(compare_ttest(:,1));node_labels(compare_ttest(:,2))]';
t_names=array2table(node_names, "VariableNames", ["Node 1 Name","Node 2 Name"]);
ttest_results_conn = [t_names table];


% Calculate FDR 
pvals=reshape(table2array(ttest_results_conn(:,7)),[],ncomp);
[~,ncomparisons]=size(pvals);
fdr=zeros(size(pvals));
for i=1:ncomparisons
    [~, q] = mafdr(pvals(:,i));
    fdr(:,i)=q;
end

fdr=reshape(fdr,[],1);
fdr_table=array2table(fdr, "VariableNames", ["q-value"]);
ttest_results_conn = [t_names table fdr_table];



end

