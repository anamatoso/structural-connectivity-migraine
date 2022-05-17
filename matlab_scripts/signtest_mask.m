function [significance_mask] = signtest_mask(a)
% a is matrix of size 116x116xn_people. Comes from load data.
% significance_mask is 1 is connection is significant and 0 otherwise

significance_mask=zeros(116,116);
for i=1:116
    for j=1:116
        data=squeeze(a(i,j,:));
        k=kstest(data);
        if k==1
            [~,h] = signtest(data,0,'Tail','right');
            
        else
            [h,~] = ttest2(data,0,'Tail','right');
        end
        if isnan(h)
            h=1;
        end
        significance_mask(i,j)=h;
    end
end
end


