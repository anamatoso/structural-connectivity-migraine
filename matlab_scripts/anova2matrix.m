function [matrix] = anova2matrix(anovacompare)

table=table2array(anovacompare(:,3:8));

for i=1:length(table)
    matrix(table(i,1),table(i,2))=table(i,6)*table(i,5);
    matrix(table(i,2),table(i,1))=table(i,6)*table(i,5);
end
end

