function [metrics_labels] = get_label_metrics()
% This funtion prints the labels of each metric
BC=strings(1,116);
Ci=strings(1,116);
EC=strings(1,116);
Qi=strings(1,116);
D=strings(1,116);
RC=strings(1,115);
node_labels=get_label_nodes("AAL116_labels.txt");
for i=1:116
    if i<=115
        BC(i)="BC_"+node_labels(i);
        Ci(i)="Ci_"+node_labels(i);
        EC(i)="EC_"+node_labels(i);
        Qi(i)="Qi_"+node_labels(i);
        D(i)="D_"+node_labels(i);
        RC(i)="RC_"+num2str(i);
        
    else
        BC(i)="BC_"+node_labels(i);
        Ci(i)="Ci_"+node_labels(i);
        EC(i)="EC_"+node_labels(i);
        Qi(i)="Qi_"+node_labels(i);
        D(i)="D_"+node_labels(i);
    end
end

metrics_labels=[BC "L" "GE" "C" Ci EC Qi "Q" RC D "mean_D" "T" "S" "A"];
end

