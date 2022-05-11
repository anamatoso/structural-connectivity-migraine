%% Load data from matrices
format long
% directory where connectivity matrices are
dir='/Users/ana/Documents/Ana/universidade/Tese/Code/matlab_scripts/matrix_data';
dir_roi='/Users/ana/Documents/Ana/universidade/Tese/Code/matlab_scripts/roi_sizes';

% Controls midcyle
HC_midcycle_mrtrix=load_data_mrtrix(dir,'*midcycle*mrtrix*bval2.csv'); %116 x 116 x n_people
%HC_midcycle_fsl=load_data_fsl(dir,'*midcycle*fsl*',dir_roi);

% Controls premenstrual
HC_premenstrual_mrtrix=load_data_mrtrix(dir,'*premenstrual*mrtrix*bval2.csv'); %116 x 116 x n_people
%HC_midcycle_fsl=load_data_fsl(dir,'*premenstrual*fsl*',dir_roi);

% Patients interictal
M_interictal_mrtrix=load_data_mrtrix(dir,'*interictal*mrtrix*bval2.csv');
%M_interictal_fsl=load_data_fsl(dir,'*interictal*fsl*',dir_roi);

% Patients ictal
M_ictal_mrtrix=load_data_mrtrix(dir,'*-ictal*mrtrix*bval2.csv');
%M_ictal_fsl=load_data_fsl(dir,'*ictal*fsl*',dir_roi);

% Note: they are not normalized by number of streamlines
connectomes={HC_midcycle_mrtrix M_interictal_mrtrix HC_premenstrual_mrtrix M_ictal_mrtrix};
n_conditions=length(connectomes);

% Calculate people per situation
n_people=zeros(1,n_conditions);
for i=1:n_conditions
    conmat=connectomes{i};
    s=size(conmat);
    n_people(i)=s(end);
end

clear dir s conmat i dir_roi
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
%load("metrics.mat")
%% Calculate metrics
% 24s per person
clear metrics
for i=1:n_conditions
    conmats=connectomes{i};
    for p=1:n_people(i)
        mat=conmats(:,:,p); % connectivity matrix
        m(:,p)=calculate_metrics(mat);
    end
    metrics{i}=m;
end
clear i p mat conmats m m2
%% Analysis of results
[n_metrics,~]=size(metrics{1});
metrics_labels=get_label_metrics();

metrics_intervsmid=[];
metrics_icvspre=[];

for m=1:n_metrics
    hc_mid=metrics{1}(m,:);
    mig_inter=metrics{2}(m,:);
    hc_pre=metrics{3}(m,:);
    mig_ic=metrics{4}(m,:);
    if (isequal(hc_mid,zeros(1,length(hc_mid))) && isequal(mig_inter,zeros(1,length(mig_inter)))) || isempty(hc_mid) || isempty(mig_inter)
        continue
    end
    x = [hc_mid mig_inter];
    g = [zeros(1,length(hc_mid)),ones(1,length(mig_inter))];
   
    
    if ttest2(hc_mid,mig_inter)==1
        metrics_intervsmid=[metrics_intervsmid m];
    end
    
    if (isequal(hc_pre,zeros(1,length(hc_pre))) && isequal(mig_ic,zeros(1,length(mig_ic)))) || isempty(hc_pre) || isempty(mig_ic)
        continue
    end
    x = [hc_pre mig_ic];
    g = [zeros(1,length(hc_pre)),ones(1,length(mig_ic))];
    if ttest2(hc_pre,mig_ic)==1
        metrics_icvspre=[metrics_icvspre m];
    end
    
    
    
end
node_labels = get_label_nodes("AAL116_labels.txt");
label_metrics_sign=metrics_labels(metrics_intervsmid);
clear hc mig x g

%% Visualization of results - Rich club

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
        case 3
            color='c';
        case 4
            color='y';
    end
    
    plot(linspace(1,115,115),mean(met(469:583,:),2),color)
    title("Rich club coefficient curve")
    xlabel("Degree")
    ylabel("Rich club coefficient")
    hold on
end

% plot significance
for m=469:583
    hc_mid=metrics{1}(m,:);
    mig_inter=metrics{2}(m,:);
    
    %     if (isequal(hc,zeros(1,length(hc))) && isequal(mig,zeros(1,length(mig)))) || isempty(hc) || isempty(mig)
    %         continue
    %     end
    
    if ttest2(hc_mid,mig_inter)==1
        plot(m-468,1,"+k")
        hold on
    end
end
hold off

clear m hc mig i met color
%% Visualization of results - Plot General metrics L, GE, C, Q, T, S,A
idx=[117 118 119 468 700 701 702];
for i=1:length(idx)
    index=idx(i);
    figure;
    
    hc_mid=metrics{1}(index,:);
    mig_inter=metrics{2}(index,:);
    
    x = [hc_mid mig_inter];
    g = [zeros(1,length(hc_mid)),ones(1,length(mig_inter))];
    
    boxplot(x,g,'Labels',{'HC','M'})
    title(metrics_labels(index))
    ylabel(metrics_labels(index))
    
    if ttest2(hc_mid,mig_inter)==1
        hold on
        text(2.1,1.01*quantile(mig_inter,0.75), '*','FontSize',15,'Color','black');
    end
    hold off
end

clear i idx index hc mig x g

%% Visualization of results - Node Strength
figure;
x=["HC M"];

title("Node Strength")
xlabel("Node")
ylabel("Strength")
n_index=[];
node=[];
bar(metrics_mean(642:688,:)); hold on
%plot significance
for m=584:699
    hc_mid=metrics{1}(m,:);
    mig_inter=metrics{2}(m,:);
    
    if (isequal(hc_mid,zeros(1,length(hc_mid))) && isequal(mig_inter,zeros(1,length(mig_inter)))) || isempty(hc_mid) || isempty(mig_inter)
        continue
    end
    
    if ttest2(hc_mid,mig_inter)==1
        text(m-583,max(mean(hc_mid),mean(mig_inter))+0.1, '*','FontSize',15); hold on
        n_index=[n_index m];
        node=[node m-583];
    end
end

%xticks([59 60 61 62 67 68 69 70 77 78 83 89 90 92 93 94 95 96 98 99 100 101 102 103 105])
hold off
% legend({'HC', 'M'})
% legend('Location','best')

%% Analysis of results - ANOVA - General metrics L, GE, C, Q, T, S,A
compare_anova=zeros(1,4);

for m=1:n_metrics
    hc_mid=metrics{1}(m,:);
    mig_inter=metrics{2}(m,:);
    hc_pre=metrics{3}(m,:);
    mig_ict=metrics{4}(m,:);
    
    x = [hc_mid mig_inter hc_pre mig_ict];
    g = [zeros(1,length(hc_mid)),ones(1,length(mig_inter)),2.*ones(1,length(hc_pre)),3.*ones(1,length(mig_ict))];
    
    if any(isnan(x))
        continue
    end
    
    [~, ~, stats] = anova1(x,g,'off');
    c=multcompare(stats,'ctype','bonferroni','display','off');
    
    for idx_p=1:length(c(:,end))
        if c(idx_p,end)<0.05
            compare_anova=[compare_anova;m c(idx_p,1) c(idx_p,2) c(idx_p,end)];
        end
    end
end
compare_anova=compare_anova(2:end,:);
table=array2table(compare_anova, "VariableNames", ["Metric index","Group 1", "Group 2", "P-value"]);
metrics_names=metrics_labels(compare_anova(:,1))';
t_names=array2table(metrics_names, "VariableNames", ["Metric Name"]);
ANOVA_results = [t_names table];
clear m hc_mid mig_inter hc_pre mig_ict x g c idx_p compare_anova table metrics_names t_names stats

