function [matrix] = anova2matrix(anovacompare,sign)

table=table2array(anovacompare(:,3:9));

if sign=="pos"
    for i=1:length(table)
        if table(i,7)>0
            matrix(table(i,1),table(i,2))=table(i,6);
            matrix(table(i,2),table(i,1))=table(i,6);
        else
            matrix(table(i,1),table(i,2))=0;
            matrix(table(i,2),table(i,1))=0;
        end
    end
    
else
    for i=1:length(table)
        if table(i,7)<0
            matrix(table(i,1),table(i,2))=-table(i,6);
            matrix(table(i,2),table(i,1))=-table(i,6);
        else
            matrix(table(i,1),table(i,2))=0;
            matrix(table(i,2),table(i,1))=0;
        end
    end
end
end

