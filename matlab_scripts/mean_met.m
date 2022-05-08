function [mean_metrics] = mean_met(metrics)
% This function calculates mean of the metrics
[~, n]=size(metrics);
mean_metrics=zeros(length(metrics{1}), n);

for i=1:length(metrics)
    condition=metrics{i};
    mean_metrics(:,i)=mean(condition,2);
end
end

