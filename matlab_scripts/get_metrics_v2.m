function [metrics] = get_metrics_v2(connectomes,v)
% Calculates metrics for a given set of connectomes

n_conditions=numel(connectomes);

% Calculate people per situation
n_people=zeros(1,n_conditions);
for i=1:n_conditions
    conmat=connectomes{i};
    s=size(conmat);
    n_people(i)=s(end);
end

metrics=cell(size(connectomes));
topprogress=sum(n_people);
progress=0;
textprogressbar('calculating metrics: ')
for i=1:n_conditions
    conmats=connectomes{i};
    clear m
    for p=1:n_people(i)
        mat=conmats(:,:,p); % connectivity matrix
        %disp(i+","+p)
        m(:,p)=calculate_metrics_v2(mat,v);
        progress=progress+1;
        textprogressbar(progress/topprogress*100);
    end
    metrics{i}=m;
end
textprogressbar('done');
end