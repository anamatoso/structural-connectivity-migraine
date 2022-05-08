%% Load data from matrices
format long
% directory where connectivity matrices are
dir='/Users/ana/Documents/Ana/universidade/Tese/Code/matlab_scripts/matrix_data';

% Controls midcyle
HC_midcycle_mrtrix=load_data(dir,'*midcycle*mrtrix*bval2.csv'); %116 x 116 x n_people
%HC_midcycle_fsl=load_data(dir,'*midcycle*fsl*');

% Patients interictal
M_interictal_mrtrix=load_data(dir,'*interictal*mrtrix*bval2.csv');
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
        mat=conmats(:,:,p); % connectivity matrix
        m(:,p)=calculate_metrics(mat);
        %m2(:,p)=calculate_metrics_v2(mat);
    end
    metrics{i}=m;
    %metrics2{i}=m2;
end
clear i p mat conmats m m2
%% Analysis of results
[n_metrics,~]=size(metrics{1});
metrics_labels=get_label_metrics();

metrics_sign=[];

for m=1:n_metrics
    hc=metrics{1}(m,:);hc=hc(~isnan(hc));
    mig=metrics{2}(m,:);mig=mig(~isnan(mig));
    
    if (isequal(hc,zeros(1,length(hc))) && isequal(mig,zeros(1,length(mig)))) || isempty(hc) || isempty(mig)
        continue
    end
    x = [hc mig];
    g = [zeros(1,length(hc)),ones(1,length(mig))];
    
    if ttest2(hc,mig)
        metrics_sign=[metrics_sign m];
        
        %figure
        %boxplot(x,g,'Labels',{'HC','M'})
        %hold on
        %text(2.1,1.01*quantile(mig,0.75), '*','FontSize',14,'Color','red');
        %hold off
    end
    
    
end
disp(metrics_labels(metrics_sign))
%% Visualization of results

metrics_mean=mean_met(metrics);

% Rich club coefficient
figure;
for i=1:n_conditions
    met=metrics{i};
    switch i
        case 1
            color='b';
        case 2
            color='r';
    end
    for p=1:n_people(i)
        plot(linspace(1,115,115),met(469:583,p),color)
        title("Rich club coefficient curve")
        xlabel("Degree")
        ylabel("Rich club coefficient")
        hold on
    end
end

% L, GE, C, Q , T, S
idx=[117 118 119 468 700 701];
for i=1:length(idx)
    index=idx(i);
    figure;
    
    hc=metrics{1}(index,:);hc=hc(~isnan(hc));
    mig=metrics{2}(index,:);mig=mig(~isnan(mig));
    
    x = [hc mig];
    g = [zeros(1,length(hc)),ones(1,length(mig))];
    
    boxplot(x,g,'Labels',{'HC','M'})
    title(metrics_labels(index))
    ylabel(metrics_labels(index))
   
    if ttest2(hc,mig)
        hold on
        text(2.1,1.01*quantile(mig,0.75), '*','FontSize',25);
        hold off
    end
end







