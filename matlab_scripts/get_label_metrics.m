function [metrics_labels] = get_label_metrics()
% This funtion prints the labels of each metric
BC=strings(1,116);
Ci=strings(1,116);
EC=strings(1,116);
Qi=strings(1,116);
D=strings(1,116);
RC=strings(1,115);
for i=1:116
    if i<=115
        BC(i)="BC_"+num2str(i);
        Ci(i)="Ci_"+num2str(i);
        EC(i)="EC_"+num2str(i);
        Qi(i)="Qi_"+num2str(i);
        D(i)="D_"+num2str(i);
        RC(i)="RC_"+num2str(i);
        
    else
        BC(i)="BC_"+num2str(i);
        Ci(i)="Ci_"+num2str(i);
        EC(i)="EC_"+num2str(i);
        Qi(i)="Qi_"+num2str(i);
        D(i)="D_"+num2str(i);
    end
end

metrics_labels=[BC "L" "GE" "C" Ci EC Qi "Q" RC D "T" "S"];
end

