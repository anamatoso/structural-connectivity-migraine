%% Load data from matrices
clear all
close all
format long
% directory where connectivity matrices are
dir='/Users/ana/Documents/Ana/universidade/Tese/Code/matlab_scripts/matrix_data';
dir_roi='/Users/ana/Documents/Ana/universidade/Tese/Code/matlab_scripts/roi_sizes';

atlas="AAL116"; 
if atlas=="AAL116" pattern="_intersect"; else pattern="*"+atlas; end

threshold=000;

% Controls midcyle
HC_midcycle_mrtrix=load_data(dir,"*midcycle*mrtrix*bval2"+pattern+".csv",dir_roi, "mrtrix",threshold); 
HC_midcycle_fsl=load_data(dir,"*midcycle*fsl*bval2",dir_roi, "fsl",threshold);

% Controls premenstrual
HC_premenstrual_mrtrix=load_data(dir,"*premenstrual*mrtrix*bval2"+pattern+".csv",dir_roi, "mrtrix",threshold);
% HC_premenstrual_fsl=load_data_fsl(dir,"*premenstrual*fsl*bval2",dir_roi, "fsl",threshold);

% Patients interictal
M_interictal_mrtrix=load_data(dir,"*interictal*mrtrix*bval2"+pattern+".csv",dir_roi, "mrtrix",threshold);
M_interictal_fsl=load_data(dir,"*interictal*fsl*bval2",dir_roi, "fsl",threshold);

% Patients ictal
M_ictal_mrtrix=load_data(dir,"*-ictal*mrtrix*bval2"+pattern+".csv",dir_roi, "mrtrix",threshold);
% M_ictal_fsl=load_data_fsl(dir,"*ictal*fsl*bval2",dir_roi, "fsl",threshold);

%connectomes={HC_midcycle_mrtrix HC_midcycle_fsl HC_premenstrual_mrtrix;M_interictal_mrtrix M_interictal_fsl M_ictal_mrtrix};
connectomes={HC_midcycle_mrtrix M_interictal_mrtrix};

patient_labels=["HC-midcycle-mrtrix" "M-interictal-mrtrix" "HC-midcycle-fsl" "M-interictal-fsl" "HC-premenstrual-mrtrix" "M-ictal-mrtrix"];
n_conditions=numel(connectomes);

% Calculate people per situation
n_people=zeros(1,n_conditions);
for i=1:n_conditions
    conmat=connectomes{i};
    s=size(conmat);
    n_people(i)=s(end);
end
node_labels=get_label_nodes(atlas+"_labels.txt");

% figure("color","w");imagesc(connectomes{3}(:,:,4));colorbar;colormap jet

clear pattern dir s conmat i dir_roi M_interictal_fsl HC_midcycle_mrtrix HC_midcycle_fsl HC_premenstrual_mrtrix M_interictal_mrtrix M_ictal_mrtrix

%% Compare matrices and matrix entries
% Plot matrices
figure('color','w');
subplot(1,2,1);imagesc(connectomes{2}(:,:,9)); colormap jet;colorbar;title("MRTrix")
subplot(1,2,2);imagesc(connectomes{4}(:,:,9)); colormap jet;colorbar;title("FSL")

% Plot boxplots
figure('color','w');
subplot(1,2,1);
boxplot([reshape(connectomes{1}(:,:,9),116*116,1) reshape(connectomes{3}(:,:,9),116*116,1)],"Labels",["MRTrix" "FSL"]);title("Distribution of matrix values")
annotation('rectangle',[0.18 0.14 0.24 0.1],'Color','green')
subplot(1,2,2);
boxplot([reshape(connectomes{1}(:,:,9),116*116,1) reshape(connectomes{3}(:,:,9),116*116,1)],"Labels",["MRTrix" "FSL"]);title("Zoom")
ylim([-0.02e-3 0.5e-3])

% Plot histograms
figure('color','w');
subplot(2,2,1);histogram(connectomes{1}(:,:,:),40);title("Distribution of matrix values in MRtrix HC");set(gca, 'YScale', 'log')
subplot(2,2,2);histogram(connectomes{2}(:,:,:),40);title("Distribution of matrix values in MRtrix M");set(gca, 'YScale', 'log')
subplot(2,2,3);histogram(connectomes{3}(:,:,:),40);title("Distribution of matrix values in FSL HC");set(gca, 'YScale', 'log')
subplot(2,2,4);histogram(connectomes{4}(:,:,:),40);title("Distribution of matrix values in FSL M");set(gca, 'YScale', 'log')

%% Compute differences between matrices
nnodes=length(node_labels);
%X = randi(nnodes,4,2);
difference_matrix=zeros(nnodes,nnodes);
difference_matrix_p=ones(nnodes,nnodes);
for i=1:nnodes-1
    for j=i+1:nnodes
            mrtrix=squeeze(connectomes{1}(i,j,:))';
            fsl=squeeze(connectomes{3}(i,j,:))';
            [h,p] = ttest2(mrtrix,fsl,"alpha",0.05/(116*115/2));
            difference_matrix(i,j)=h; difference_matrix(j,i)=h;
            difference_matrix_p(i,j)=p; difference_matrix_p(j,i)=p;
    end
end



connectomes2=connectomes;
for i=1:n_conditions
    for p=1:n_people(i)
        connectomes2{i}(:,:,p)=connectomes{i}(:,:,p).*(ones(size(difference_matrix(:,:)))-difference_matrix(:,:));
    end
end

figure('color','w');
subplot(1,2,1);imagesc(difference_matrix); colormap jet;colorbar;title("Coherence between fsl and mrtrix")
subplot(1,2,2);imagesc(difference_matrix_p); colormap jet;colorbar;title("Coherence between fsl and mrtrix-pvalue")


figure('color','w','Position', [100 100 2000 1000]);
subplot(2,2,1);imagesc(connectomes{1}(:,:,9)); colormap jet;colorbar;title("MRTrix")
subplot(2,2,2);imagesc(connectomes{3}(:,:,9)); colormap jet;colorbar;title("FSL")
subplot(2,2,3);imagesc(connectomes2{1}(:,:,9)); colormap jet;colorbar;title("MRTrix only coherent")
subplot(2,2,4);imagesc(connectomes2{3}(:,:,9)); colormap jet;colorbar;title("FSL only coherent")
clear mrtrix fsl x g X difference_matrix i j h p

%% Scatter plot and correlation coefficients
nnodes=length(node_labels);
scatterv=[[],[]];
for c=1:2
    for p=1:n_people(c)
        for i=1:nnodes-1
            for j=i+1:nnodes
                scatterv=[scatterv; connectomes{c}(i,j,p) connectomes{c+2}(i,j,p)];
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

mdl = fitlm(scatterv(:,1),scatterv(:,2));
disp("r^2 with normal scale: " +num2str(mdl.Rsquared.Ordinary))
mdllog = fitlm(scattervlog(:,1),scattervlog(:,2));
disp("r^2 with log-log scale: " +num2str(mdllog.Rsquared.Ordinary))
figure('color','w');
subplot(1,2,1);plot(mdl);text(6.5e-3,0.25e-3,"R^2="+num2str(mdl.Rsquared.Ordinary))
title("Normal Scale");xlabel("Mrtrix");ylabel("FSL");
subplot(1,2,2);plot(mdllog);text(-7.5,-2.5,"R^2="+num2str(mdllog.Rsquared.Ordinary))
title("Logarithm Scale");xlabel("Mrtrix");ylabel("FSL");

disp("Normal Scale")
R= corr(scatterv,"Type","Pearson");
disp("Pearson's correlation coefficient: "+ num2str(R(1,2)))
R= corr(scatterv,"Type","Kendall");
disp("Kendall's correlation coefficient: "+ num2str(R(1,2)))
R= corr(scatterv,"Type","Spearman");
disp("Spearman's correlation coefficient: "+ num2str(R(1,2)))
disp(" ")
disp("Logarithmic Scale")
R= corr(scattervlog,"Type","Pearson");
disp("Pearson's correlation coefficient: "+ num2str(R(1,2)))
R= corr(scattervlog,"Type","Kendall");
disp("Kendall's correlation coefficient: "+ num2str(R(1,2)))
R= corr(scattervlog,"Type","Spearman");
disp("Spearman's correlation coefficient: "+ num2str(R(1,2)))

clear R npeople nnodes scatterv c p i j scattervlog todelete mdl mdllog R

%% Analyse only a subnetwork of the connectome (Optional)
subnetwork=[1:90];

for i=1:n_conditions
    connectomes{i}=connectomes{i}(subnetwork,subnetwork,:);
end

node_labels=node_labels(subnetwork);
clear i subnetwork

%% Remap matrix (Optional)

% aggregate adjacent zones
idx_map=[1 2 repmat([3 4],1,7) 5 6 7 8 9 10 repmat([3 4],1,2) 11 12 13 14 repmat([15 16],1,3)...
    17 18 17 18 19 20 21 22 23 24 25 26 repmat([27 28],1,3) 29 30 1 2 repmat([31 32],1,2)...
    33 34 35 36 37 38 1 2 39 40 41 42 43 44 45 46 47 48 repmat([49 50],1,5) 51*ones(1,4) 52*ones(1,14) 53*ones(1,8)];
node_labels=["precentralL" "precentralR" "FrontalL" "FrontalR" "RolandicL" "RolandicR" "suppmotorL" "suppmotorR"...
    "olfactoryL" "olfactoryR" "rectusL" "rectusR" "insulaL" "insulaR" "CingulumL" "cingulumR"...
    "hippocampusL" "hippocampusR" "amygdalaL" "amygdalaR" "calcarineL" "calcarineR" "cuneusL" "cuneusR"...
    "LingualL" "LingualR" "occipitalL" "occipitalR" "fusiforeL" "fusiformR" "parietalL" "parietalR"...
    "suppramarginalL" "suppramarginalR" "angularL" "angularR" "precuneusL" "precuneusR" "caudateL" "caudateR"...
    "putamenL" "putamenR" "pallidiumL" "pallidiumR" "thalamusL" "thalamusR" ...
    "heschlL" "heschlR" "temporalL" "temporalR" "Cerebelum_crus" "Cerebelum_crus" "vermis"];
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
[n_nodes,~,~]=size(connectomes{1});
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

version_metrics=3;%  1=nodal metrics, 2=general metrics
clear metrics
metrics=cell(size(connectomes));
for i=1:n_conditions
    conmats=connectomes{i};
    clear m
    for p=1:n_people(i)
        mat=conmats(:,:,p); % connectivity matrix
        disp(i+","+p)
        m(:,p)=calculate_metrics(mat,version_metrics);
    end
    metrics{i}=m;
end

clear i p mat conmats m m2

%% Analysis of results
version_metrics=3;
metrics_labels=get_label_metrics(version_metrics,node_labels);
comparisons=[1 2];

ttest_results = ttest_compare_v2(metrics,metrics_labels,version_metrics,length(node_labels),comparisons);

%writetable(ttest_results, 'ttest_results.xlsx');
clear comparisons

%% Visualization of results: metrics
idx_metrics=[1:7];
idx_groups=[1:2];
patient_labels=["HC-midcycle-mrtrix" "M-mrtrix" "HC-fsl" "M-fsl" "HC-premenstrual-mrtrix" "M-ictal-mrtrix"];

plot_boxplots(metrics,idx_metrics,idx_groups,metrics_labels,patient_labels,version_metrics,116)

clear idx_metrics idx_groups
%% For visualization in BrainNet nodes AAL116
% nodal metrics:
nodestrength=(349:464); bc=(1:116); lC=(117:232); ec=(233:348);
m=[nodestrength;bc;lC;ec];
names=["nodestrength" "bc" "lC" "ec"];

for i=1:4
    qvalues=table2array(ttest_results(m(i,:),8)); 
    diff=table2array(ttest_results(m(i,:),7));
    nodes_degree_color = nodes_color_size(qvalues,diff,0.05,node_labels);
    nodefile = table(makenodefile("aal116_MNIcoord.txt",node_labels,nodes_degree_color));
    writetable(nodefile, names(i)+'_significantq_fsl_sym.txt','Delimiter',' ','WriteVariableNames', 0);
end
clear pvalues diff nodes_degree_color nodefile nodestrength bc lC ec i
%% Analysis of connectivity
comparisons=[1 2];
ttest_results_conn = ttest_compare_conn(connectomes,node_labels,comparisons);

%writetable(ANOVA_results_conn, 'ttest_results_conn.xlsx');
clear comparisons
%% Visualization of results connectivity
idx_metrics=[61 63;1 11;1 35;1 52];
idx_groups=[3 4];
plot_boxplots_conn(connectomes,idx_metrics,idx_groups,patient_labels,node_labels)

%% For visualization in BrainNet edges AAL116
matrix = anova2matrix(ttest_results_conn,"neg");
writematrix(matrix, 'edges_AAL116_tscore_neg.txt','Delimiter',' ');

%% Determine hub nodes
mean_matrices=calculate_mean_matrix(connectomes);
%node_labels=node_labels(subnetwork);
n_nodes=length(mean_matrices);
hubnodes=strings(round(n_nodes*0.2),4); % SO 2 CONDITIONS
for i=1:2 % SO 2 CONDITIONS
    [idx, values]=hub_nodes(mean_matrices(:,:,i)); %indices of nodes
    labels=node_labels(idx); %names of nodes
    hubnodes(:,2*(i)-1:2*i)=[labels' values'];
end
hubnodestable=array2table(hubnodes,'VariableNames',{'HC_midcycle','HC_midcycle_BC','M_interictal','M_interictal_BC'});%,'HC_premenstrual','M_ictal'});

color_size = hubnodes_color_size(hubnodes(:,1),hubnodes(:,2));

matrix = makenodefile(color_size);
T=table(matrix);
writetable(T, 'hubnodes_HC.txt','Delimiter',' ','WriteVariableNames', 0);
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


