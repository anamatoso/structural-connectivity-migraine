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

[allconnectomes,~,~,node_labels,~] = get_data(dir,dir_roi,atlas,threshold,1:4,false);

%load("allconnectomes.mat")

clear threshold atlas dir dir_roi

%% Transform to N4 vs N2

N2N4connectomes={allconnectomes{2} allconnectomes{4}};% cell 1x2, in which each cell is a 1x2 cell in which each cell is each group

for g=1:2
    connectomes=N2N4connectomes{g};
    connectomes={connectomes{1},connectomes{2}};%,connectomes{5},connectomes{6}};
    N2N4connectomes{g}=connectomes;
end

clear connectomes g allconnectomes

%% Calculate metrics
version=2;
allmetrics_norms=cell(1,2);
allmetrics_groups=cell(1,2);
for norm=1:2
    allmetrics_groups{norm}=get_metrics_v2(N2N4connectomes{norm},version);
    allmetrics_norms{norm}=cat(2,allmetrics_groups{norm}{1},allmetrics_groups{norm}{2});
end

metrics_labels=get_label_metrics(version,node_labels);

clear norm version

%% Plot results comparing both normalizations

%Create array for p values
table1=zeros(length(metrics_labels),2);

% Iterate through metrics
for metric=1:length(metrics_labels)

    % Calculate p values
    N2=allmetrics_norms{1}(metric,:);
    N4=allmetrics_norms{2}(metric,:);
    [~,p]=ttest2(N2,N4,'vartype', 'unequal'); % not tottally correct but stated
    table1(metric,:)=[metric p];
    
    % Plot figure
    figure('color','w')
    boxplot([N2' N4'],'Labels',["N1","N2"])
    title(metrics_labels(metric),'interpreter', 'none','FontSize',20,'FontWeight','normal','FontName','Arial')
    text(2.1,1*max(max(N2),max(N4)), "p="+num2str(p),'FontSize',12,'Color','black');
    set(gca,'FontSize',20)
end

% Create final table and display the metrics where p<0.05
T1=array2table(metrics_labels','VariableNames',"Metric Name");
T2=array2table(table1,'VariableNames',["Metric Index","P-value"]);
T=[T1 T2];
disp(T(T.("P-value")<0.05,:))

clear N2 N4 metric p T1 T2

%% See significance between each group comparison in both norms 

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
T=[T2(:,1) T1 T2(:,2:4)];
disp(T(T.("P-value")<0.05,:))

clear metric p_1 p_2 table T1 T2 

%% Plot bloxchart for each metric
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

clear metric xlabels g group norm xgroup x colourdata i datapoint positionaldata n_people metric_data

%% Analysis of results

comparisons=[1 2];

ttest_results=cell(size(allmetrics_groups));
for g=1:length(allmetrics_groups)
    metrics=allmetrics_groups{g};
    ttest_results{g} = ttest_compare_v2(metrics,metrics_labels,version,length(node_labels),comparisons);
end
%writetable(ttest_results, 'ttest_results.xlsx');
clear comparisons g metrics

%% For visualization in BrainNet nodes AAL116 - Groups

norm=1;
bc=(1:116); lC=(117:232); ec=(233:348); nodestrength=(349:464);
m=[bc;lC;ec;nodestrength];
names=["BC" "Ci" "EC" "D"];
ttest_results2=ttest_results{norm};
for g=1:4
    pvalues=table2array(ttest_results2(m(g,:),5)); 
    diff=table2array(ttest_results2(m(g,:),7));
    nodes_degree_color = nodes_color_size(pvalues,diff,0.05);
    nodefile = table(makenodefile("aal116_MNIcoord.txt",node_labels,nodes_degree_color));
    writetable(nodefile, 'nodes/ismrm23/'+names(g)+'_midinter_n'+string(norm)+'.txt','Delimiter',' ','WriteVariableNames', 0);
end

clear diff nodes_degree_color nodefile nodestrength bc lC ec g pvalues m diff norm names

%% For visualization in BrainNet nodes AAL116 - Normalisations
names=["BC" "Ci" "EC" "D"];
node_labels=get_label_nodes("AAL116_labels_number.txt");
for g=1:4
    pvalues=table1(m(g,:),2); 
    diff=ones(116,1); % all changed (no matter positively or negatively)
    nodes_degree_color = nodes_color_size(pvalues,diff,0.05/116/4);
    nodefile = table(makenodefile("aal116_MNIcoord.txt",node_labels,nodes_degree_color));
    writetable(nodefile, 'nodes/ismrm23/'+names(g)+'_diffnorms.txt','Delimiter',' ','WriteVariableNames', 0);
end

clear diff nodes_degree_color nodefile nodestrength bc lC ec g pvalues m diff names

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

clear roi_size idx fullFileName baseFileName k theFiles filePattern F dir_roi mid int
