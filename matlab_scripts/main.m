%% Load data from matrices
clear variables
close all
format long
% directory where connectivity matrices are
dir='/Users/ana/Documents/Ana/universidade/Tese/Code/matlab_scripts/matrix_data';
dir_roi='/Users/ana/Documents/Ana/universidade/Tese/Code/matlab_scripts/roi_sizes';

atlas="AAL116"; 
if atlas=="AAL116"; pattern="_intersect"; else; pattern="*"+atlas; end

threshold=000;
for norm=1:4
    % Controls midcyle
    HC_midcycle_mrtrix=load_data(dir,"*midcycle*mrtrix*bval2"+pattern+".csv",dir_roi, "mrtrix",threshold,norm); 
    HC_midcycle_fsl=load_data(dir,"*midcycle*fsl*bval2_omat3",dir_roi, "fsl",threshold,norm);
    
    % Controls premenstrual
    HC_premenstrual_mrtrix=load_data(dir,"*premenstrual*mrtrix*bval2"+pattern+".csv",dir_roi, "mrtrix",threshold,norm);
    HC_premenstrual_fsl=load_data(dir,"*premenstrual*fsl*bval2_omat3",dir_roi, "fsl",threshold,norm);
    
    % Patients interictal
    M_interictal_mrtrix=load_data(dir,"*interictal*mrtrix*bval2"+pattern+".csv",dir_roi, "mrtrix",threshold,norm);
    M_interictal_fsl=load_data(dir,"*interictal*fsl*bval2_omat3",dir_roi, "fsl",threshold,norm);
    
    % Patients ictal
    M_ictal_mrtrix=load_data(dir,"*-ictal*mrtrix*bval2"+pattern+".csv",dir_roi, "mrtrix",threshold,norm);
    M_ictal_fsl=load_data(dir,"*-ictal*fsl*bval2_omat3",dir_roi, "fsl",threshold,norm);

    allconnectomes{norm}={HC_midcycle_mrtrix HC_midcycle_fsl HC_premenstrual_mrtrix HC_premenstrual_fsl;...
    M_interictal_mrtrix M_interictal_fsl M_ictal_mrtrix M_ictal_fsl};
end
%connectomes={HC_midcycle_mrtrix M_interictal_mrtrix};
%connectomes=allconnectomes{1};
n_conditions=numel(allconnectomes{1});

% Calculate people per situation
n_people=zeros(1,n_conditions);
for i=1:n_conditions
    conmat=allconnectomes{1}{i};
    s=size(conmat);
    n_people(i)=s(end);
end
node_labels=get_label_nodes(atlas+"_labels.txt");
condition_names=["MRtrix HC midcycle" "MRtrix M interictal" "FSL HC midcycle" "FSL M interictal" "MRtrix M premenstrual" "MRtrix M ictal" "FSL M premenstrual" "FSL M ictal"];

% figure("color","w");imagesc(connectomes{3}(:,:,4));colorbar;colormap jet

clear pattern dir s norm conmat i dir_roi M_interictal_fsl HC_midcycle_mrtrix HC_midcycle_fsl HC_premenstrual_mrtrix HC_premenstrual_fsl M_interictal_mrtrix M_ictal_mrtrix M_ictal_fsl

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
clear x1 x2

% Plot boxplots
% figure('color','w');
% subplot(1,2,1);
% boxplot([reshape(connectomes{1}(:,:,9),116*116,1) reshape(connectomes{3}(:,:,9),116*116,1)],"Labels",["MRTrix" "FSL"]);title("Distribution of matrix values")
% annotation('rectangle',[0.18 0.14 0.24 0.1],'Color','green')
% subplot(1,2,2);
% boxplot([reshape(connectomes{1}(:,:,9),116*116,1) reshape(connectomes{3}(:,:,9),116*116,1)],"Labels",["MRTrix" "FSL"]);title("Zoom")
% ylim([-0.02e-3 0.5e-3])

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
nnodes=length(node_labels);
connectomes=allconnectomes{1};
scatterv=[[],[]];
topprogress=(sum(n_people)/2)*((nnodes*nnodes-nnodes)/2);
progress=0;
textprogressbar('Getting values: ')
for c=[1 2 5 6]%linspace(1,7,4)
    for p=1:n_people(c)
        for i=1:nnodes-1
            for j=i+1:nnodes
                scatterv=[scatterv; connectomes{c}(i,j,p) connectomes{c+2}(i,j,p)];
                progress=progress+1;
                textprogressbar(progress/topprogress*100);
            end
        end
    end
end

scattervlog=log10(scatterv);
todelete=[];
for i=1:length(scattervlog)
    if isinf(scattervlog(i,1)) || isinf(scattervlog(i,2))
        todelete=[todelete i];
    end
end
todelete=flip(todelete);
for i=1:length(todelete)
    scattervlog(todelete(i),:)=[];
end

figure('color','w')
set(gca,'FontSize',20)
scatterhist(scattervlog(:,1),scattervlog(:,2), 'NBins',[40,40],'Direction','out',Marker='x'); xlabel('MRTrix');ylabel('FSL')
%%
mdl = fitlm(scatterv(:,1),scatterv(:,2));
disp("r^2 with normal scale: " +num2str(mdl.Rsquared.Ordinary))
mdllog = fitlm(scattervlog(:,1),scattervlog(:,2));
disp("r^2 with log-log scale: " +num2str(mdllog.Rsquared.Ordinary))

figure('color','w');
subplot(1,2,1);plot(mdl);text(6.5e-3,0.25e-3,"R^2="+num2str(mdl.Rsquared.Ordinary))
title("Linear Scale");xlabel("Mrtrix");ylabel("FSL");
subplot(1,2,2);plot(mdllog);text(-7.5,-2.5,"R^2="+num2str(mdllog.Rsquared.Ordinary))
title("Logarithmic Scale");xlabel("Mrtrix");ylabel("FSL");

disp("Normal Scale")
R= corr(scatterv,"Type","Pearson");
disp("Pearson's correlation coefficient: "+ num2str(R(1,2)))
% R= corr(scatterv,"Type","Kendall");
% disp("Kendall's correlation coefficient: "+ num2str(R(1,2)))
R= corr(scatterv,"Type","Spearman");
disp("Spearman's correlation coefficient: "+ num2str(R(1,2)))
disp(" ")
disp("Logarithmic Scale")
R= corr(scattervlog,"Type","Pearson");
disp("Pearson's correlation coefficient: "+ num2str(R(1,2)))
% R= corr(scattervlog,"Type","Kendall");
% disp("Kendall's correlation coefficient: "+ num2str(R(1,2)))
R= corr(scattervlog,"Type","Spearman");
disp("Spearman's correlation coefficient: "+ num2str(R(1,2)))

clear R npeople nnodes c p i j mdl mdllog R

%% Analyse only a subnetwork of the connectome (Optional)
subnetwork=1:90;

for i=1:n_conditions
    connectomes{i}=connectomes{i}(subnetwork,subnetwork,:);
end

node_labels=node_labels(subnetwork);
clear i subnetwork

%% Remap matrix (Optional)
connectomes=allconnectomes{1};
% aggregate adjacent zones
idx_map=[1 2 repmat([3 4],1,7) 5 6 7 8 9 10 repmat([3 4],1,2) 11 12 13 14 repmat([15 16],1,3)...
    17 18 17 18 19 20 21 22 23 24 25 26 repmat([27 28],1,3) 29 30 1 2 repmat([31 32],1,2)...
    33 34 35 36 37 38 1 2 39 40 41 42 43 44 45 46 47 48 repmat([49 50],1,5) repmat([51 52],1,2) repmat([53 54],1,7) 55*ones(1,8)];
node_labels=["precentralL" "precentralR" "FrontalL" "FrontalR" "RolandicL" "RolandicR" "suppmotorL" "suppmotorR"...
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
connectomes=newconnectome;
clear i p newconnectome idx_map
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
%connectomes=rescale_connectomes(connectomes,n_people);
% connectomes =connectome2aal90(connectomes);

version_metrics=2;%  1=nodal metrics, 2=general metrics
%load("allmetrics"+version_metrics+".mat")

allmetrics=cell(size(allconnectomes));
for i=1:length(allconnectomes)
    connectomes=allconnectomes{i};
    allmetrics{i}=get_metrics(connectomes,version_metrics);
end

clear i
%% Analysis of results
metrics_labels=get_label_metrics(version_metrics,node_labels);
comparison_HCvsP=[1 2;3 4;5 6;7 8];
comparison_MRtrixvsFSL=[1 3;2 4;5 7;6 8];
comparison_cycle=[1 5;3 7;2 6;4 8];

comparisons=comparison_HCvsP;

ttest_results=cell(size(allmetrics));
for i=1:length(allmetrics)
    metrics=allmetrics{i};
    ttest_results{i} = ttest_compare_v2(metrics,metrics_labels,version_metrics,length(node_labels),comparisons);
end
%writetable(ttest_results, 'ttest_results.xlsx');
clear comparisons comparison_HCvsP comparison_MRtrixvsFSL comparison_cycle i

%% Visualization of results: metrics
idx_metrics=[1];
idx_groups=[3 4 7 8];
metrics=allmetrics{1};
condition_names=["MRtrix-HC-midcycle" "MRtrix-M-interictal" "FSL-HC-midcycle" "FSL-M-interictal" "MRtrix-M-premenstrual" "MRtrix-M-ictal" "FSL-M-premenstrual" "FSL-M-ictal"];

%plot_boxplots(metrics,idx_metrics,idx_groups,metrics_labels,condition_names,version_metrics,116)
plot_boxplots_both(metrics,idx_metrics,[],0)
clear idx_metrics idx_groups
%% For visualization in BrainNet nodes AAL116
% nodal metrics:
nodestrength=(349:464); bc=(1:116); lC=(117:232); ec=(233:348);
m=[nodestrength;bc;lC;ec];
names=["D" "BC" "Ci" "EC"];
ttest_results2=ttest_results{1};
for i=1:4
    qvalues=table2array(ttest_results2(m(i,:),5)); 
    diff=table2array(ttest_results2(m(i,:),7));
    nodes_degree_color = nodes_color_size(qvalues,diff,0.05);
    nodefile = table(makenodefile("aal116_MNIcoord.txt",node_labels,nodes_degree_color));
    writetable(nodefile, 'nodes/cycle/prob/'+names(i)+'_FSL_intervsict.txt','Delimiter',' ','WriteVariableNames', 0);
end
clear pvalues diff nodes_degree_color nodefile nodestrength bc lC ec i qvalues pvals m diff
%% Analysis of connectivity
comparison_HCvsP=[1 2;3 4;5 6;7 8];
comparison_MRtrixvsFSL=[1 3;2 4;5 7;6 8];
comparison_cycle=[1 5;3 7;2 6;4 8];
connectomes=allconnectomes{1};

comparisons=[4 8];
ttest_results_conn = ttest_compare_conn(connectomes,node_labels,comparisons);

%writetable(ANOVA_results_conn, 'ttest_results_conn.xlsx');
clear comparisons comparison_HCvsP comparison_MRtrixvsFSL comparison_cycle
%% Visualization of results connectivity
idx_metrics=[61 63;1 11;1 35;1 52];
idx_groups=[3 4];
plot_boxplots_conn(connectomes,idx_metrics,idx_groups,patient_labels,node_labels)

%% Create matrix to plot
n_edges=25;

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

clear i j idx mat P1 P2

%% For visualization in BrainNet edges AAL116
matrix = anova2matrix(ttest_results_conn,"neg");
writematrix(matrix, 'edges_AAL116_tscore_neg.txt','Delimiter',' ');

%% Determine hub nodes
mean_matrices=calculate_mean_matrix(connectomes);
%node_labels=node_labels(subnetwork);
n_nodes=length(mean_matrices);
hubnodes=cell(2,4); % SO 2 CONDITIONS
for i=1:8 % SO 2 CONDITIONS
    [idx, values]=hub_nodes(mean_matrices(:,:,i)); %indices of nodes
    labels=node_labels(idx); %names of nodes
    hubnodes{i}=[labels' values'];
end
%hubnodestable=array2table(hubnodes,'VariableNames',{'HC_midcycle','HC_midcycle_BC','M_interictal','M_interictal_BC'});%,'HC_premenstrual','M_ictal'});

%color_size = hubnodes_color_size(hubnodes(:,1),hubnodes(:,2));

%matrix = makenodefile(color_size);
%T=table(matrix);
%writetable(T, 'hubnodes_HC.txt','Delimiter',' ','WriteVariableNames', 0);
clear i idx labels T matrix n_nodes values

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

%% Visualization of results - Node Strength
figure;
x=["HC M"];
metrics_mean=mean_met(metrics);

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


