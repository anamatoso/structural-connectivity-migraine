function [metrics_labels] = get_label_metrics(a,node_labels)
% This funtion prints the labels of each metric given a choice of metrics and the label of each node

n_nodes=length(node_labels);
if a==1
    BC=strings(1,n_nodes);
    Ci=strings(1,n_nodes);
    EC=strings(1,n_nodes);
    D=strings(1,n_nodes);
    RC=strings(1,n_nodes-1);
    
    for i=1:n_nodes
        if i<=n_nodes-1
            BC(i)="BC_"+node_labels(i);
            Ci(i)="Ci_"+node_labels(i);
            EC(i)="EC_"+node_labels(i);
            D(i)="D_"+node_labels(i);
            RC(i)="RC_"+num2str(i);
            
        else
            BC(i)="BC_"+node_labels(i);
            Ci(i)="Ci_"+node_labels(i);
            EC(i)="EC_"+node_labels(i);
            D(i)="D_"+node_labels(i);
        end
    end
    metrics_labels=[BC Ci EC RC D];
    
elseif a==2
    metrics_labels=["L" "GE" "C" "Q" "mean_D" "T" "S"];
  
end

