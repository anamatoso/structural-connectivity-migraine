function [x,y] = extract_groups(data,groups,idx1,idx2)


x=zeros(1,sum(groups(:) == idx1));
y=zeros(1,sum(groups(:) == idx2));

ix=1;iy=1;
for i=1:length(data)
    if groups(i)==idx1
        x(ix)=data(i);
        ix=ix+1;
    elseif groups(i)==idx2
        y(iy)=data(i);
        iy=iy+1;
    end 

end
end