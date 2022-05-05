%% Load data from matrices
% directory where connectivity matrices are
dir='/Users/ana/Documents/Ana/universidade/Tese/Code/matlab_scripts/matrix_data';

% Controls midcyle
HC_midcycle_mrtrix=load_data(dir,'*midcycle*mrtrix*'); %116 x 116 x n_people
%HC_midcycle_fsl=load_data(dir,'*midcycle*fsl*'); 

% Patients interictal
M_interictal_mrtrix=load_data(dir,'*interictal*mrtrix*');
%M_interictal_fsl=load_data(dir,'*interictal*fsl*');

% Note: they are not normalized by number of streamlines
connectomes={HC_midcycle_mrtrix M_interictal_mrtrix}; 
n_conditions=length(connectomes);

% Calculate people per situation
n_people=zeros(1,n_conditions);
for i=1:n_conditions
    conmat=connectomes{i};
    s=size(conmat);
    n_people(i)=s(end);
end

clear dir s conmat i
%% Test for spurious connections

significance_mask=zeros(116,116,n_conditions);
for i=1:n_conditions
    significance_mask(:,:,i) = signtest_mask(connectomes{i});
    for p=1:n_people(i)
        connectomes{i}(:,:,p)=connectomes{i}(:,:,p).*significance_mask(:,:,i);
    end
end
% imagesc(significance_mask(:,:,i)); colormap jet
clear i p
%% Calculate metrics- load result of script below
load("metrics.mat")
%% Calculate metrics
% 24s per person
for i=1:n_conditions
    conmats=connectomes{i};
    for p=1:n_people(i)
        mat=conmats(:,:,p);                 % connectivity matrix
        %m(:,p)=calculate_metrics(mat);
        m2(:,p)=calculate_metrics_v2(mat);
    end
    %metrics{i}=m;
    metrics2{i}=m2;
end
clear i p mat conmats m m2
%% Analysis of results
[n_metrics,~]=size(metrics2{1});
metrics_labels=get_label_metrics();

metrics_sign=[];
for m=1:n_metrics
    hc=metrics2{1}(m,:);hc=hc(~isnan(hc));
    mig=metrics2{2}(m,:);mig=mig(~isnan(mig));
    if (isequal(hc,zeros(1,length(hc))) && isequal(mig,zeros(1,length(mig)))) || isempty(hc) || isempty(mig)
        continue
    end
    x = [hc mig];
    g = [zeros(1,n_people(1)),ones(1,n_people(2))];
    
     if ttest2(hc,mig, 'Alpha', 0.05/4) || ttest2(hc,mig, 'Alpha', 0.05/4,'Tail','right') || ttest2(hc,mig, 'Alpha', 0.05/4,'Tail','left')
%         figure
%         boxplot(x,g,'Labels',{'HC','M'})
%         hold on
%         text(2.1,1.01*quantile(mig,0.75), '*','FontSize',14,'Color','red');
%         hold off
        metrics_sign=[metrics_sign m];
    end
    
end
disp(metrics_labels(metrics_sign))
