%% Compare Normalizations/Algorithms/conditions
clc
clear variables
close all
format long
%% Load data from matrices
dir=strcat(pwd,'/matrix_data_prob');
dir_roi=strcat(pwd,'/roi_sizes');
atlas="AAL116";
threshold=0;

[allconnectomes,~,~,node_labels,~] = get_data(dir,dir_roi,atlas,threshold,1:4);

%load("allconnectomes.mat")

clear dir dir_roi atlas threshold

%% Select only N2 and N4

N2N4connectomes={allconnectomes{2} allconnectomes{4}};% cell 1x2, in which each cell is a 1x2 cell in which each cell is each group

for g=1:2
    connectomes=N2N4connectomes{g};
    connectomes={connectomes{1},connectomes{2}};%,connectomes{5},connectomes{6}};
    N2N4connectomes{g}=connectomes;
end

clear allconnectomes g connectomes

%% Calculate metrics

version=2; % Choose type of metrics to calculate

allmetrics_norms=cell(1,2);
allmetrics_groups=cell(1,2);
for norm=1:2
    allmetrics_groups{norm}=get_metrics_v2(N2N4connectomes{norm},version);
    allmetrics_norms{norm}=cat(2,allmetrics_groups{norm}{1},allmetrics_groups{norm}{2});
end

metrics_labels=get_label_metrics(version,node_labels);

clear norm version

%% See results comparing both normalizations and plot

show=1; % Plot boxplots

%Create array for p values
table=zeros(length(metrics_labels),2);

% Iterate through metrics
for metric=1:length(metrics_labels)

    % Calculate p values
    N2=allmetrics_norms{1}(metric,:);
    N4=allmetrics_norms{2}(metric,:);
    [~,p]=ttest2(N2,N4); % not tottally correct but stated
    table(metric,:)=[metric p];
    
    if show
        % Plot figure
        figure('color','w')
        boxplot([N2' N4'],'Labels',["N1","N2"])
        title(metrics_labels(metric),'interpreter', 'none','FontSize',20,'FontWeight','normal','FontName','Arial')
        text(2.1,1*max(max(N2),max(N4)), "p="+num2str(p),'FontSize',12,'Color','black');
        set(gca,'FontSize',20)
    end
end

% Create final table and display the metrics where p<0.05
T1=array2table(metrics_labels','VariableNames',"Metric Name");
T2=array2table(table,'VariableNames',["Metric Index","P-value"]);
T_norms=[T1 T2];
disp(T_norms(T_norms.("P-value")<0.05,:))

clear N2 N4 metric p T1 T2 table show

%% For visualization in BrainNet nodes AAL116 - Normalisations

% Define metrics
bc=(1:116); lC=(117:232); ec=(233:348); nodestrength=(349:464);
m=[bc;lC;ec;nodestrength];
names=["BC" "Ci" "EC" "D"];

% Iterate through metrics
for metric=1:4
    pvalues=table2array(T_norms(m(metric,:),3)); 
    diff=ones(116,1); % all changed (no matter positively or negatively)
    nodes_degree_color = nodes_color_size(pvalues,diff,0.05/116/2);
    nodefile = table(makenodefile("aal116_MNIcoord.txt",node_labels,nodes_degree_color));
    writetable(nodefile, 'nodes/ismrm23/'+names(metric)+'_diffnorms.txt','Delimiter',' ','WriteVariableNames', 0);
end

clear bc lC ec nodestrength m names metric pvalues diff nodes_degree_color nodefile

%% See results comparing each group in both normalisations 

%Create array for p values
table=zeros(2*length(metrics_labels),4);

% Iterate through metrics
for metric=1:length(metrics_labels)

    % Calculate p values
    [p_1,~]=ranksum(allmetrics_groups{1}{1}(metric,:),allmetrics_groups{1}{2}(metric,:));
    [p_2,~]=ranksum(allmetrics_groups{2}{1}(metric,:),allmetrics_groups{2}{2}(metric,:));
    
    % Insert in table
    table(2*metric-1,:)  =[1 metric p_1 median(allmetrics_groups{1}{2}(metric,:))-median(allmetrics_groups{1}{1}(metric,:))];
    table(2*metric,:)=[2 metric p_2 median(allmetrics_groups{2}{2}(metric,:))-median(allmetrics_groups{2}{1}(metric,:))];
end

% Create final table and display the metrics where p<0.05
T1=array2table(metrics_labels(table(:,2))',"VariableNames","Metric Name");
T2=array2table(table,'VariableNames',["Norm","Metric Index","P-value","Diff"]);
T_groups=[T2(:,1) T1 T2(:,2:4)];
disp(T_groups(T_groups.("P-value")<0.05,:))

clear table metric p_1 p_2 T1 T2 

%% Plot boxchart for each metric

[~,n_people]= size(allmetrics_norms{1});
for metric=1:length(metrics_labels)
    
    % Prepare data
    x=zeros(1,2*n_people);
    xgroup=zeros(1,2*n_people);
    colourdata=zeros(1,2*n_people);
    datapoint=1;
    for group=1:2
        for norm=1:2
            metric_data=allmetrics_groups{norm}{group}(metric,:);
            for i=1:length(metric_data)
                x(datapoint)=metric_data(i);
                xgroup(datapoint)=group;
                colourdata(datapoint)=norm;
                datapoint=datapoint+1;
            end
        end
    end

    % Turn xgroup into labels aka categorical variable
    xlabels={'Controls' 'Migraineurs'};
    positionaldata=discretize(xgroup,1:3,'categorical',xlabels);
    
    % Plot figure
    figure('color','w','Position',[200 100 895 500])
    boxchart(positionaldata, x,"GroupByColor",colourdata)
    grid on
    title(metrics_labels(metric),'interpreter', 'none','FontSize',20,'FontWeight','normal','FontName','Arial')
    xline(1.5)
    legend(["N1", "N2"])
    set(gca,'FontSize',20)

end

clear n_people metric x xgroup colourdata datapoint group norm metric_data i xlabels positionaldata  

%% Analysis of results - Groups

ttest_results=cell(size(allmetrics_groups));
for g=1:length(allmetrics_groups)
    metrics=allmetrics_groups{g};
    ttest_results{g} = ttest_compare_v2(metrics,metrics_labels,version,length(node_labels),1:2);
end

clear g metrics

%% For visualization in BrainNet nodes AAL116 - Groups

% Define metrics and chose normalisation
norm=2;
bc=(1:116); lC=(117:232); ec=(233:348); nodestrength=(349:464);
m=[bc;lC;ec;nodestrength];
names=["BC" "Ci" "EC" "D"];

% Iterate through metrics
ttest_results2=ttest_results{norm};
for metric=1:4
    pvalues=table2array(ttest_results2(m(metric,:),5)); 
    diff=table2array(ttest_results2(m(metric,:),7));
    nodes_degree_color = nodes_color_size(pvalues,diff,0.05);
    nodefile = table(makenodefile("aal116_MNIcoord.txt",node_labels,nodes_degree_color));
    writetable(nodefile, 'nodes/ismrm23/'+names(metric)+'_midinter_n'+string(norm)+'.txt','Delimiter',' ','WriteVariableNames', 0);
end

clear norm bc lC ec nodestrength m names metric pvalues diff nodes_degree_color nodefile

%% Plot histogram ROI sizes

dir_roi=strcat(pwd,'/roi_sizes');
F = dir_roi;
filePattern = fullfile(F,"*_intersect*");
theFiles = dir(filePattern);

roisizes={zeros(116,15), zeros(116,14)};
mid=1;int=1;
for k =1:length(theFiles)
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(theFiles(k).folder, baseFileName);
    roi_size=importdata(fullFileName); 
    if contains(fullFileName,"midcycle")
        roisizes{1}(:,mid)=roi_size;
        mid=mid+1;
    elseif contains(fullFileName,"interictal")
        roisizes{2}(:,int)=roi_size;
        int=int+1;
    end
end

% Plot bar plot for each group with average size of each region
for idx=1:2
    figure('color','w','Position',[200 100 895 500])
    bar(mean(roisizes{idx}, 2))
    xticks([10 30 50 70 90 110])
    set(gca,'FontSize',15)
    title("Group "+ idx)
    xlabel("Region index in AAL116")
    ylabel("Average Region Volume (in voxels)")
end

% Plot bar plot for both groups with average size of each region
figure('color','w','Position',[200 100 895 500])
bar(mean(cat(2,roisizes{1},roisizes{2}),2))
xticks([10 30 50 70 90 110])
set(gca,'FontSize',15)
title("Both Groups")
xlabel("Region index in AAL116")
ylabel("Average Region Volume (in voxels)")

clear dir_roi F filePattern theFiles mid int k baseFileName fullFileName roi_size roi_size
