%% Load data from matrices
clear all
format long
% directory where connectivity matrices are
dir='/Users/ana/Documents/Ana/universidade/Tese/Code/matlab_scripts/matrix_data';
dir_roi='/Users/ana/Documents/Ana/universidade/Tese/Code/matlab_scripts/roi_sizes';

atlas="AAL116"; 
if atlas=="AAL116" pattern="_intersect"; else pattern="*"+atlas; end

% Controls midcyle
HC_midcycle_mrtrix=load_data(dir,"*midcycle*mrtrix*bval2"+pattern+".csv",dir_roi, "mrtrix"); 
HC_midcycle_fsl=load_data(dir,"*midcycle*fsl*",dir_roi, "fsl");

% Controls premenstrual
HC_premenstrual_mrtrix=load_data(dir,"*premenstrual*mrtrix*bval2"+pattern+".csv",dir_roi, "mrtrix");
% HC_premenstrual_fsl=load_data_fsl(dir,"*premenstrual*fsl*",dir_roi, "fsl");

% Patients interictal
M_interictal_mrtrix=load_data(dir,"*interictal*mrtrix*bval2"+pattern+".csv",dir_roi, "mrtrix");
M_interictal_fsl=load_data(dir,"*interictal*fsl*",dir_roi, "fsl");

% Patients ictal
M_ictal_mrtrix=load_data(dir,"*-ictal*mrtrix*bval2"+pattern+".csv",dir_roi, "mrtrix");
% M_ictal_fsl=load_data_fsl(dir,"*ictal*fsl*",dir_roi, "fsl");

connectomes={HC_midcycle_mrtrix HC_midcycle_fsl HC_premenstrual_mrtrix;M_interictal_mrtrix M_interictal_fsl M_ictal_mrtrix};
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
clear pattern dir s conmat i dir_roi M_interictal_fsl HC_midcycle_mrtrix HC_midcycle_fsl HC_premenstrual_mrtrix M_interictal_mrtrix M_ictal_mrtrix
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
%% Remove spurious connections
%connectomes=rescale_connectomes(connectomes,n_people);
[n_nodes,~,~]=size(connectomes{1});
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
%connectomes=rescale_connectomes(connectomes,n_people);
% connectomes =connectome2aal90(connectomes);

version_metrics=3;%  1=all metrics, 2=general metrics
clear metrics
metrics=cell(size(connectomes));
for i=1:n_conditions
    conmats=connectomes{i};
    for p=1:n_people(i)
        mat=conmats(:,:,p); % connectivity matrix
        m(:,p)=calculate_metrics(mat,version_metrics);
    end
    metrics{i}=m;
    clear m
end


clear i p mat conmats m m2

%% Analysis of results
version_metrics=3;
metrics_labels=get_label_metrics(version_metrics,node_labels);
comparisons=[1 2];
ttest_results = ttest_compare(metrics,metrics_labels,version_metrics,length(node_labels),comparisons);

%writetable(ttest_results, 'ANOVA_results.xlsx');
clear comparisons

%% Visualization of results: metrics
idx_metrics=[103 63 331 325];
idx_groups=[1 2];

plot_boxplots(metrics,idx_metrics,idx_groups,metrics_labels,patient_labels)

clear idx_metrics idx_groups
%% For visualization in BrainNet nodes AAL116
% nodal metrics:
nodestrength=(464:579); bc=(1:116);lC=(117:232);ec=(233:348);

pvalues=table2array(ttest_results(bc,5)); 
diff=table2array(ttest_results(bc,6));
nodes_degree_color = nodes_color_size(pvalues,diff,0.05,node_labels);
nodefile = table(makenodefile("aal116_MNIcoord.txt",node_labels,nodes_degree_color));
writetable(nodefile, 'bc_significant.txt','Delimiter',' ','WriteVariableNames', 0);
clear pvalues diff nodes_degree_color nodefile nodestrength bc lC ec
%% Analysis of connectivity
comparisons=[1 2];
ttest_results_conn = ttest_compare_conn(connectomes,node_labels,comparisons);

%writetable(ANOVA_results_conn, 'ttest_results_conn.xlsx');
clear comparisons
%% Visualization of results connectivity
idx_metrics=[27 64;5 74;28 54;84 95];
idx_groups=[1 2];
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


