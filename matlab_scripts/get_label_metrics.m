function [metrics_labels] = get_label_metrics(a,node_labels)
% This funtion prints the labels of each metric given a choice of metrics and the label of each node

n_nodes=length(node_labels);
if a==1
    RC=strings(1,n_nodes-1);
    
    for i=1:n_nodes-1     
        RC(i)="RC_"+num2str(i);
    end
    
    metrics_labels=RC;
    
elseif a==2
    metrics_labels=["L" "GE" "C" "Q" "mean_D" "T" "S"];
    
elseif a==3
    BC=strings(1,n_nodes);
    Ci=strings(1,n_nodes);
    EC=strings(1,n_nodes);
    D=strings(1,n_nodes);
    
    for i=1:n_nodes
        BC(i)="BC_"+node_labels(i);
        Ci(i)="Ci_"+node_labels(i);
        EC(i)="EC_"+node_labels(i);
        D(i)="D_"+node_labels(i);
        
    end
    metrics_labels=[BC Ci EC D];
    
    
end

