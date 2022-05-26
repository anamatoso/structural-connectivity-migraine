function [ANOVA_results_conn] = anova_compare_conn(connectomes,node_labels,varargin)

compare_anova=zeros(1,6);
n_nodes=length(node_labels);
n_comparisons=(n_nodes*n_nodes-n_nodes)/2;

if ~isempty(varargin)
    if varargin{1}=="True"
        threshold=inf;
    else
        threshold=0.05/n_comparisons;
    end
end
for m=1:n_nodes-1
    for n=m+1:n_nodes
        hc_mid=squeeze(connectomes{1}(m,n,:))';
        mig_inter=squeeze(connectomes{2}(m,n,:))';
        %hc_pre=squeeze(connectomes{3}(m,n,:))';
        %mig_ict=squeeze(connectomes{4}(m,n,:))';
        
        x = [hc_mid mig_inter];% hc_pre mig_ict];
        g = [zeros(1,length(hc_mid)),ones(1,length(mig_inter))];%,2.*ones(1,length(hc_pre)),3.*ones(1,length(mig_ict))];
        
        if any(isnan(x)) || all(x==0)
            continue
        end
        
        [~, ~, stats] = anova1(x,g,'off');
        c=multcompare(stats,'display','off');
        
        for idx_p=1:length(c(:,end))
            
            if c(idx_p,end)<threshold  
                switch c(idx_p,4)>0
                    case 1
                        compare_anova=[compare_anova;m n c(idx_p,1) c(idx_p,2) c(idx_p,end)*n_comparisons 1];
                    case 0
                        compare_anova=[compare_anova;m n c(idx_p,1) c(idx_p,2) c(idx_p,end)*n_comparisons -1];
                end
            end
        end
    end
end
compare_anova=compare_anova(2:end,:);
table=array2table(compare_anova, "VariableNames", ["Node 1","Node 2","Group 1", "Group 2", "P-value","Difference"]);
node_names=[node_labels(compare_anova(:,1));node_labels(compare_anova(:,2))]';
t_names=array2table(node_names, "VariableNames", ["Node 1 Name","Node 2 Name"]);
ANOVA_results_conn = [t_names table];
end

