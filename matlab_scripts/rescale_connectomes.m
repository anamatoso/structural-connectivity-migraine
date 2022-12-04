function [c] = rescale_connectomes(connectomes)
% This function resclales all connectomes to [0,1]

% Determine people per group
npeople=zeros(size(connectomes));
for i=1:numel(connectomes)
    [~,~,n]=size(connectomes{i});
    npeople(i)=n;
end

% Rescale every connectome
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

