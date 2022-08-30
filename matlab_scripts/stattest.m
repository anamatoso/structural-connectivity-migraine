function [p] = stattest(data1,data2)


if ~kstest(data1) && ~kstest(data2) && ~vartest2(data1,data2)   % normal with equal variances
    [~,p] = ttest2(data1,data2);
elseif ~kstest(data1) && ~kstest(data2) && vartest2(data1,data2)% normal with unequal variances
    [~,p] = ttest2(data1,data2,'vartype', 'unequal');
else                                                            % not normal with unequal variances
    p = ranksum(data1,data2);
end
end