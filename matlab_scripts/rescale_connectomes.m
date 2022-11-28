function [c] = rescale_connectomes(connectomes,npeople)
% This function resclales all connectomes to [0,1]

c=cell(size(connectomes));
for i=1:length(connectomes)
    matrices=connectomes{i};
    newmat=zeros(size(matrices));
    for p=1:npeople(i)
        newmat(:,:,p)=rescale(matrices(:,:,p),0,1);
    end
    c{i}=newmat;
end

end

