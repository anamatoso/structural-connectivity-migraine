%% Load data from matrices
clear all
format long
% directory where connectivity matrices are
dir='/Users/ana/Documents/Ana/universidade/Tese/Code/matlab_scripts/matrix_data';
dir_roi='/Users/ana/Documents/Ana/universidade/Tese/Code/matlab_scripts/roi_sizes';

% Controls midcyle
HC_midcycle_mrtrix=load_data_mrtrix(dir,'*midcycle*mrtrix*bval2.csv'); %116 x 116 x n_people
%HC_midcycle_fsl=load_data_fsl(dir,'*midcycle*fsl*',dir_roi);
%HC_midcycle_mrtrix=cat(3,HC_midcycle_mrtrix,HC_midcycle_mrtrix);

% Controls premenstrual
HC_premenstrual_mrtrix=load_data_mrtrix(dir,'*premenstrual*mrtrix*bval2.csv'); %116 x 116 x n_people
%HC_midcycle_fsl=load_data_fsl(dir,'*premenstrual*fsl*',dir_roi);
%HC_premenstrual_mrtrix=cat(3,HC_premenstrual_mrtrix,HC_premenstrual_mrtrix);

% Patients interictal
M_interictal_mrtrix=load_data_mrtrix(dir,'*interictal*mrtrix*bval2.csv');
%M_interictal_fsl=load_data_fsl(dir,'*interictal*fsl*',dir_roi);
%M_interictal_mrtrix=cat(3,M_interictal_mrtrix,M_interictal_mrtrix);

% Patients ictal
M_ictal_mrtrix=load_data_mrtrix(dir,'*-ictal*mrtrix*bval2.csv');
%M_ictal_fsl=load_data_fsl(dir,'*ictal*fsl*',dir_roi);
%M_ictal_mrtrix=cat(3,M_ictal_mrtrix,M_ictal_mrtrix);

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
node_labels=get_label_nodes("AAL116_labels.txt");
clear dir s conmat i dir_roi HC_midcycle_mrtrix HC_midcycle_fsl HC_premenstrual_mrtrix M_interictal_mrtrix M_ictal_mrtrix
%% Analyse only a subnetwork of the connectome (Optional)
subnetwork=[1 2 3 4 9 10 29 30 31 32 57 58 59 60 69 70 71 72 73 74 77 78];
subnetwork=[1:90];
for i=1:n_conditions
    connectomes{i}=connectomes{i}(subnetwork,subnetwork,:);
end
[n_nodes,~,~]=size(connectomes{1});
node_labels=get_label_nodes("AAL116_labels.txt");
node_labels=node_labels(subnetwork);
clear i

%% Remap matrix (Optional)

% agregate L-R and all vermis 
idx_map=[];
for i=1:54
 idx_map=[idx_map i*ones(1,2)];  
end
idx_map=[idx_map 55*ones(1,8)];

% idx_map=[ones(1,2) 2*ones(1,4) 3*ones(1,4) 4*ones(1,6) 5 5 6 6 7 7 8*ones(1,4) 9 9 10 10 ...
%     11*ones(1,6) 12*ones(1,4) 13 13 14 14 15 15 16 16 17*ones(1,6) 18 18 19 19 20*ones(1,4) ...
%     21 21 22 22 23 23 24 24 25 25 26 26 27 27 28 28 29 29 30*ones(1,10) ...
%     31*ones(1,4) 32*ones(1,14) 33*ones(1,8)];
newconnectome=cell(size(connectomes));
for i=1:n_conditions
    for p=1:n_people(i)
        newconnectome{i}(:,:,p)=remap_matrix(connectomes{i}(:,:,p),idx_map);
    end
end
connectomes=newconnectome;
clear i p newconnectome idx_map
%% Remove spurious connections
%connectomes=rescale_connectomes(connectomes,n_people);
significance_mask=zeros(n_nodes,n_nodes,n_conditions);
for i=1:n_conditions
    significance_mask(:,:,i) = signtest_mask(connectomes{i});
    for p=1:n_people(i)
        connectomes{i}(:,:,p)=connectomes{i}(:,:,p).*significance_mask(:,:,i);
    end
end
% imagesc(significance_mask(:,:,1)); colormap jet;colorbar
clear i p

%% Calculate metrics
% connectomes=rescale_connectomes(connectomes,n_people);
% connectomes =connectome2aal90(connectomes);

version_metrics=2;%  1=703 metrics, 2=124 metrics, 3=8 metrics
clear metrics
for i=1:n_conditions
    conmats=connectomes{i};
    for p=1:n_people(i)
        mat=conmats(:,:,p); % connectivity matrix
        m(:,p)=calculate_metrics(mat,version_metrics);
    end
    metrics{i}=m;
end


clear i p mat conmats m m2 version_metrics

%% Set variables to the result analysis

version_metrics=2;
[n_metrics,~]=size(metrics{1});
metrics_labels=get_label_metrics(version_metrics,node_labels);

%% Analysis of results - ANOVA

ANOVA_results = anova_compare(metrics,metrics_labels,1);

%% Analysis of connectivity between nodes- ANOVA
ANOVA_results_conn = anova_compare_conn(connectomes,node_labels,1);

%% Determine hub nodes
mean_matrices=calculate_mean_matrix(connectomes);
node_labels = get_label_nodes("AAL116_labels.txt");
%node_labels=node_labels(subnetwork);
hubnodes=strings(round(n_nodes*0.2),2); % SO 2 CONDITIONS
for i=1:2 % SO 2 CONDITIONS
    idx=hub_nodes(mean_matrices(:,:,i)); %indices of nodes
    labels=node_labels(idx); %names of nodes
    hubnodes(:,i)=labels';
end
hubnodestable=array2table(hubnodes,'VariableNames',{'HC_midcycle','M_interictal'});%,'HC_premenstrual','M_ictal'});

color_size = hubnodes_color_size(hubnodes(:,1),hubnodes(:,2));

matrix = makenodefile(color_size);
T=table(matrix);
writetable(T, 'hubnodes_HC.txt','Delimiter',' ','WriteVariableNames', 0);
clear i idx labels T matrix

%% Test alternative to anova
compare_anova=zeros(1,5);

for m=1:n_metrics
    hc_mid=metrics{1}(m,:);
    mig_inter=metrics{2}(m,:);
    
    x = [hc_mid mig_inter];
    %g = [zeros(1,length(hc_mid)),ones(1,length(mig_inter))];
    
    if any(isnan(x))
        continue
    end
    
    % check variance
    switch vartest2(hc_mid,mig_inter)
        case 1 % different varience
            [h,p]=ttest2(hc_mid,mig_inter,'alpha',0.05/n_metrics,"vartype",'unequal');
        case 0 % equal varience
            [h,p]=ttest2(hc_mid,mig_inter,'alpha',0.05/n_metrics);
    end
    
    
    if h==1 % if significant difference
        c=mean(hc_mid)-mean(mig_inter);
        switch c>0 % difference is positive?
            case 1
                compare_anova=[compare_anova;m 1 2 p*n_metrics 1];
            case 0
                compare_anova=[compare_anova;m 1 2 p*n_metrics -1];
        end
    end
    
end
compare_anova=compare_anova(2:end,:);
table=array2table(compare_anova, "VariableNames", ["Metric index","Group 1", "Group 2", "P-value corrected", "Difference"]);
metrics_names=metrics_labels(compare_anova(:,1))';
t_names=array2table(metrics_names, "VariableNames", ["Metric Name"]);
ANOVA_results = [t_names table];
clear m h p hc_mid mig_inter hc_pre mig_ict x g c idx_p compare_anova table metrics_names t_names stats compare_anova table metrics_names t_names

%% Visualization of results - Rich club

%metrics_mean=mean_met(metrics);

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
            color='g';
    end
    
    plot(linspace(1,115,115),nanmean(met(469:583,:),2),color)
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
legend(["HC midcycle", "M interictal", "HC premenstrual", "M ictal"])
clear m hc mig i met color
%% Visualization of results
idx1=[1 2 3 4 121 122 123 124]; % General metrics
idx2=[94 82 109 100 87]; % metrics with smaller p
plot_boxplots(metrics,idx1,metrics_labels)

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


