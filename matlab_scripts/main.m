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
% imagesc(significance_mask(i)); colormap jet
clear i p
%% Calculate metrics
% 24s per person
for i=1:n_conditions
    conmats=connectomes{i};
    for p=1:n_people(i)
        mat=conmats(:,:,p);                 % connectivity matrix
        m(:,p)=calculate_metrics(mat);     
    end
    metrics{i}=m;
end
clear i p mat conmats
%% Analysis of results
n_metrics=length(metrics{1});
metrics2analyse=[117 118 119 583 584 585];
for i=1:length(metrics2analyse)
    m=metrics2analyse(i); %index of metric
    figure
    x = [metrics{1}(m,:) metrics{2}(m,:)];
    g = [zeros(1,n_people(1)),ones(1,n_people(2))];
    boxplot(x,g,'Labels',{'HC','M'})
     hold on
    
%     if % difference is significative
%         text(2.1,1.01*quantile(metrics{2}(m,:),0.75), '*','FontSize',14,'Color','red');
%     end
    hold off
end




%%