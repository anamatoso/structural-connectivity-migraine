%% Load data from matrices
clear variables
close all
format long

%% Load data from matrices
dir=strcat(pwd,'/matrix_data_prob');
dir_roi=strcat(pwd,'/roi_sizes');

atlas="AAL116"; 
threshold=0;

[allconnectomes,n_conditions,n_people,node_labels,condition_names] = get_data(dir,dir_roi,atlas,threshold,1:4);

% figure("color","w");imagesc(connectomes{3}(:,:,4));colorbar;colormap jet

clear threshold atlas dir_roi dir

%% Compare matrices and matrix entries

% Plot matrices
figure('color','w');
x1=1;
x2=5;
for i=1:8
    if rem(i,2) %odd
        subplot(2,4,x1);imagesc(connectomes{i}(:,:,8)); colormap jet;colorbar;title(condition_names(i))
        x1=x1+1;
    else
        subplot(2,4,x2);imagesc(connectomes{i}(:,:,8)); colormap jet;colorbar;title(condition_names(i))
        x2=x2+1;
    end
end

% Plot histograms
x1=1;
x2=5;
figure('color','w');
for i=1:8
    if rem(i,2) %odd
        subplot(2,4,x1);histogram(log10(connectomes{i}(:,:,:)),40);title(condition_names(i));set(gca, 'YScale', 'log');xlim([-12 -2])
        x1=x1+1;
    else
        subplot(2,4,x2);histogram(log10(connectomes{i}(:,:,:)),40);title(condition_names(i));set(gca, 'YScale', 'log');xlim([-12 -2])
        x2=x2+1;
    end
end
clear x1 x2

%% Scatter plot and correlation coefficients

connectomes=allconnectomes{1}; % Normalization 1

nnodes=length(node_labels);
topprogress=(sum(n_people)/2)*((nnodes*nnodes-nnodes)/2);
scatterv=zeros(topprogress,2);
progress=0;
textprogressbar('Getting values: ')

% Create vector with values in linear scale
for c=[1 2 5 6] % mrtrix
    for p=1:n_people(c)
        for i=1:nnodes-1
            for j=i+1:nnodes
                scatterv(progress+1,:)=[connectomes{c}(i,j,p) connectomes{c+2}(i,j,p)];
                progress=progress+1;
                textprogressbar(progress/topprogress*100);
            end
        end
    end
end

% Create vector with values in logarithmic scale
scattervlog=log10(scatterv);
for i=length(scattervlog):-1:1
    if isinf(scattervlog(i,1)) || isinf(scattervlog(i,2))
        scattervlog(i,:)=[];
    end
end

% Plot
figure('color','w')
scatterhist(scattervlog(:,1),scattervlog(:,2), 'NBins',[40,40],'Direction','out',Marker='x'); 
xlabel('MRTrix');ylabel('FSL')
title("Normalization 1 - Logarithmic Scale")
set(gca,'FontSize',15)

clear topprogress progress nnodes connectomes c p i j

%% Fit Regression Model

% Fit Regression
mdl = fitlm(scatterv(:,1),scatterv(:,2));
disp("r^2 with normal scale: " +num2str(mdl.Rsquared.Ordinary))
mdllog = fitlm(scattervlog(:,1),scattervlog(:,2));
disp("r^2 with log-log scale: " +num2str(mdllog.Rsquared.Ordinary))

% Plot
figure('color','w');
subplot(1,2,1);plot(mdl);text(6.5e-3,0.25e-3,"R^2="+num2str(mdl.Rsquared.Ordinary))
title("Linear Scale");xlabel("Mrtrix");ylabel("FSL");
subplot(1,2,2);plot(mdllog);text(-7.5,-2.5,"R^2="+num2str(mdllog.Rsquared.Ordinary))
title("Logarithmic Scale");xlabel("Mrtrix");ylabel("FSL");

% Display correlation coefficients
disp("Normal Scale")
R= corr(scatterv,"Type","Pearson");
disp("Pearson's correlation coefficient: "+ num2str(R(1,2)))
R= corr(scatterv,"Type","Spearman");
disp("Spearman's correlation coefficient: "+ num2str(R(1,2)))
disp(" ")

disp("Logarithmic Scale")
R= corr(scattervlog,"Type","Pearson");
disp("Pearson's correlation coefficient: "+ num2str(R(1,2)))
R= corr(scattervlog,"Type","Spearman");
disp("Spearman's correlation coefficient: "+ num2str(R(1,2)))

clear R mdl mdllog scattervlog scatterv

%% Analyse only a subnetwork of the connectome (Optional)

subnetwork=1:90; % AAL90
new_connectomes=cell(size(connectomes));
for i=1:n_conditions
    new_connectomes{i}=connectomes{i}(subnetwork,subnetwork,:);
end

%connectomes=new_connectomes;
%node_labels=node_labels(subnetwork);

clear i subnetwork

%% Remap matrix (Optional)
connectomes=allconnectomes{1}; % Normalization 1

% Eg: aggregate adjacent zones
idx_map=[1 2 repmat([3 4],1,7) 5 6 7 8 9 10 repmat([3 4],1,2) 11 12 13 14 repmat([15 16],1,3)...
    17 18 17 18 19 20 21 22 23 24 25 26 repmat([27 28],1,3) 29 30 1 2 repmat([31 32],1,2)...
    33 34 35 36 37 38 1 2 39 40 41 42 43 44 45 46 47 48 repmat([49 50],1,5) repmat([51 52],1,2) repmat([53 54],1,7) 55*ones(1,8)];
new_node_labels=["precentralL" "precentralR" "FrontalL" "FrontalR" "RolandicL" "RolandicR" "suppmotorL" "suppmotorR"...
    "olfactoryL" "olfactoryR" "rectusL" "rectusR" "insulaL" "insulaR" "CingulumL" "cingulumR"...
    "hippocampusL" "hippocampusR" "amygdalaL" "amygdalaR" "calcarineL" "calcarineR" "cuneusL" "cuneusR"...
    "LingualL" "LingualR" "occipitalL" "occipitalR" "fusiformL" "fusiformR" "parietalL" "parietalR"...
    "suppramarginalL" "suppramarginalR" "angularL" "angularR" "precuneusL" "precuneusR" "caudateL" "caudateR"...
    "putamenL" "putamenR" "pallidiumL" "pallidiumR" "thalamusL" "thalamusR" ...
    "heschlL" "heschlR" "temporalL" "temporalR" "CerebelumcrusR" "CerebelumcrusL" "CerebelumL" "CerebelumR" "vermis"];
newconnectome=cell(size(connectomes));
for i=1:n_conditions
    for p=1:n_people(i)
        newconnectome{i}(:,:,p)=remap_matrix(connectomes{i}(:,:,p),idx_map);
    end
end

clear i p idx_map
%% Remove spurious connections based on sparsity
%connectomes=rescale_connectomes(connectomes,n_people);
percentile_mask=cell(size(connectomes));

for i=1:n_conditions
    percentile_mask{i} = create_percentile_mask(connectomes{i},25);
    for p=1:n_people(i)
        connectomes{i}(:,:,p)=connectomes{i}(:,:,p).*percentile_mask{i}(:,:,p);
    end
end

%imagesc(percentile_mask(:,:,2)); colormap jet;colorbar
clear i p

%% Calculate metrics

% Optional:
% connectomes=rescale_connectomes(connectomes);

version_metrics=2;%  1=rich club coefficient, 2=global metrics, 3=nodal metrics

% Calculate metrics
allmetrics=cell(size(allconnectomes));
for i=1:length(allconnectomes)
    connectomes=allconnectomes{i};
    allmetrics{i}=get_metrics(connectomes,version_metrics);
end
metrics_labels=get_label_metrics(version_metrics,node_labels);

clear i connectomes
%% Statistical Analysis of Results

comparison_HCvsP=[1 2;3 4;5 6;7 8];
%comparison_MRtrixvsFSL=[1 3;2 4;5 7;6 8];
%comparison_cycle=[1 5;3 7;2 6;4 8];

comparisons=comparison_HCvsP;

ttest_results=cell(size(allmetrics));
for i=1:length(allmetrics)
    metrics=allmetrics{i};
    ttest_results{i} = ttest_compare_v2(metrics,metrics_labels,numel(node_labels),comparisons);
end
%writetable(ttest_results, 'ttest_results.xlsx');

clear comparisons comparison_HCvsP comparison_MRtrixvsFSL comparison_cycle i metrics

%% Visualization of results: metrics
idx_metrics=1:7;
%idx_groups=[3 4 7 8];
metrics=allmetrics{1};
%condition_names=["MRtrix-HC-midcycle" "MRtrix-M-interictal" "FSL-HC-midcycle" "FSL-M-interictal" "MRtrix-M-premenstrual" "MRtrix-M-ictal" "FSL-M-premenstrual" "FSL-M-ictal"];

%plot_boxplots(metrics,idx_metrics,idx_groups,metrics_labels,condition_names)
plot_boxplots_both(metrics,idx_metrics,[],0)

clear idx_metrics idx_groups metrics condition_names
%% For visualization in BrainNet nodes AAL116
% nodal metrics:
nodestrength=(349:464); bc=(1:116); lC=(117:232); ec=(233:348);
m=[nodestrength;bc;lC;ec];
names=["D" "BC" "Ci" "EC"];

ttest_results2=ttest_results{1}; % normalization 1
for i=1:4 % for all nodal metrics
    qvalues=table2array(ttest_results2(m(i,:),5)); 
    diff=table2array(ttest_results2(m(i,:),7));
    nodes_degree_color = nodes_color_size(qvalues,diff,0.05);
    nodefile = table(makenodefile("aal116_MNIcoord.txt",node_labels,nodes_degree_color));
    writetable(nodefile, 'nodes/cycle/prob/'+names(i)+'_FSL_intervsict.txt','Delimiter',' ','WriteVariableNames', 0);
end

clear diff nodes_degree_color nodefile nodestrength bc lC ec i qvalues m ttest_results2
%% Analysis of connectivity
comparison_HCvsP=[1 2;3 4;5 6;7 8];
%comparison_MRtrixvsFSL=[1 3;2 4;5 7;6 8];
%comparison_cycle=[1 5;3 7;2 6;4 8];
connectomes=allconnectomes{1};

comparisons=comparison_HCvsP;

ttest_results_conn = ttest_compare_conn(connectomes,node_labels,comparisons);

%writetable(ANOVA_results_conn, 'ttest_results_conn.xlsx');

clear comparisons comparison_HCvsP comparison_MRtrixvsFSL comparison_cycle connectomes
%% Visualization of results connectivity
idx_metrics=[61 63;1 11;1 35;1 52];
idx_groups=[3 4];

plot_boxplots_conn(connectomes,idx_metrics,idx_groups,patient_labels,node_labels)

%% Create matrix to plot connectivity using mne-connectivity in python

n_edges=25; % Significant edges to plot

mat=table2array(ttest_results_conn);
nnodes=length(node_labels);
matrix=zeros(nnodes,nnodes);
idx=1;
for i=1:nnodes-1
    for j=i+1:nnodes
        a=(1-str2double(mat(idx,7)))*str2double(mat(idx,9));
        if abs(a)<0.95
            a=0;
        end
        if a<0
            a=a+1;
        end
        matrix(i,j)=a;
        matrix(j,i)=a;
        idx=idx+1;
    end
end

matreshape=reshape(matrix,1,[]);
matreshape=nonzeros(matreshape);
matreshape=sort(matreshape);
matreshape(2:2:end, :) = [];   
P1 = length(matreshape)-n_edges;
matreshape(P1)
%matreshape=flip(matreshape);
P2 = n_edges;
matreshape(P2)

for i=1:nnodes-1
    for j=i+1:nnodes
        if matrix(i,j)<matreshape(length(matreshape)-n_edges) && matrix(i,j)>matreshape(n_edges)
            matrix(i,j)=0;
            matrix(j,i)=0;
        end

    end
end
matrix=10*matrix;

for i=1:nnodes-1
    for j=i+1:nnodes
        if matrix(i,j)>9.5
            matrix(i,j)=matrix(i,j)-9;
            matrix(j,i)=matrix(i,j);
        end
    end
end

min(min(nonzeros(matrix)))

clear i j idx mat P1 P2 n_edges nnodes a

