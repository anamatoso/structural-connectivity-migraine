function [x,y] = extract_groups(data,groups,idx1,idx2)
% Given a data vector and a groups vector that contains the group of each
% element in data (length(data)=length(groups)), this function extracts the
% data points where groups==idx1 to x and the data points where
% groups==idx2 to y.

x=zeros(1,sum(groups(:) == idx1));
y=zeros(1,sum(groups(:) == idx2));

i_x=1;i_y=1;
for idx=1:length(data)
    if groups(idx)==idx1
        x(i_x)=data(idx);
        i_x=i_x+1;
    elseif groups(idx)==idx2
        y(i_y)=data(idx);
        i_y=i_y+1;
    end 
end
end