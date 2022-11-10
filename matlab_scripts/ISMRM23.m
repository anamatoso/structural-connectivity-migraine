%% Compare Normalizations/Algorithms/conditions
clc
clear variables
close all
format long
%% Load data from matrices

dir='/Users/ana/Documents/Ana/universidade/Tese/Code/matlab_scripts/matrix_data_prob';
dir_roi='/Users/ana/Documents/Ana/universidade/Tese/Code/matlab_scripts/roi_sizes';
atlas="AAL116";
threshold=0;
normalizations=[1 2 3 4];
[allconnectomes,n_conditions,n_people,node_labels,condition_names] = get_data(dir,dir_roi,atlas,threshold,normalizations,false);

%load("allconnectomes.mat")

clear threshold atlas dir dir_roi normalizations

%% Transform to N4 vs N2

N2N4connectomes={allconnectomes{2} allconnectomes{4}};
N2N4connectomes_v2={allconnectomes{2} allconnectomes{4}};

condition_names=[condition_names(1) condition_names(2) condition_names(5) condition_names(6)];

for i=1:2
    connectomes=N2N4connectomes{i};
    connectomes={connectomes{1},connectomes{2}};%,connectomes{5},connectomes{6}};
    connectomes=cat(3,connectomes{1},connectomes{2});%,connectomes{3},connectomes{4});
    N2N4connectomes{i}=connectomes;

    connectomes=N2N4connectomes_v2{i};
    connectomes={connectomes{1},connectomes{2}};%,connectomes{5},connectomes{6}};
    N2N4connectomes_v2{i}=connectomes;
end

clear connectomes i
%% Calculate metrics
version=2;
allmetrics=get_metrics_v2(N2N4connectomes,version);
allmetrics_v2=cell(1,2);
for norm=1:2
    allmetrics_v2{norm}=get_metrics_v2(N2N4connectomes_v2{norm},version);
end

metrics_labels=get_label_metrics(version,node_labels);

clear norm

%% Plot results comparing all groups
table1=zeros(116,4);
for metric=1:length(metrics_labels)
    N2=allmetrics{1}(metric,:);
    N4=allmetrics{2}(metric,:);
    [~,p]=ttest2(N2,N4);
    table1(metric)=p;
    %disp(metrics_labels(metric)+": "+p)
    
    figure('color','w')
    set(gca,'FontSize',15)
    boxplot([N2' N4'],'Labels',["N1","N2"])
    title(metrics_labels(metric),'interpreter', 'none','FontSize',20,'FontWeight','normal','FontName','Arial')
    text(2.1,1*max(max(N2),max(N4)), "p="+num2str(p),'FontSize',12,'Color','black');
    set(gca,'FontSize',20)
end
clear N2 N4 metric p
%% See significance between each group comparison in both norms 
people=[15 14 15 9];
comparison=[1 2];
table=zeros(2*4*length(metrics_labels),5);
i=1;
for metric=1:length(metrics_labels)
    %disp(metrics_labels(metric))
    for comp=1:1
        groups=comparison(comp,:); %[1 3]
        %disp("   "+condition_names(groups(1))+" vs "+condition_names(groups(2)))
        [p_1,~]=ranksum(allmetrics_v2{1}{groups(1)}(metric,:),allmetrics_v2{1}{groups(2)}(metric,:));
        [p_2,~]=ranksum(allmetrics_v2{2}{groups(1)}(metric,:),allmetrics_v2{2}{groups(2)}(metric,:));
        
        table(i,:)=[1 metric 10*groups(1)+groups(2) p_1 median(allmetrics_v2{1}{groups(2)}(metric,:))-median(allmetrics_v2{1}{groups(1)}(metric,:))];
        table(i+1,:)=[2 metric 10*groups(1)+groups(2) p_2 median(allmetrics_v2{2}{groups(2)}(metric,:))-median(allmetrics_v2{2}{groups(1)}(metric,:))];
        i=i+2;
        if p_1<0.05
            disp(metrics_labels(metric)+", "+condition_names(groups(1))+" vs "+condition_names(groups(2))+", N2: "+p_1)
        elseif p_2<0.05
            disp(metrics_labels(metric)+", "+condition_names(groups(1))+" vs "+condition_names(groups(2))+", N4: "+p_2)
        end
    end
end
T=array2table(table,'VariableNames',["Norm","Metric","Comparison","P-value","Diff"]);

clear metric comp p_1 p_2 groups i table
%% plot bloxchart for each metric

for metric=1:7
    figure('color','w','Position',[386 100 895 500])
    

    x=[];
    xgroup=[];
    i=1;
    for group=1:2
        for norm=1:2
            x=[x allmetrics_v2{norm}{group}(metric,:)];
            xgroup=[xgroup i*ones(1,people(group))];
        end
        i=i+1;
    end

    xlabels={'Controls' 'Migraineurs'};% 'HC-Premenstrual' 'M-Ictal'};
    positionaldata=discretize(xgroup,1:3,'categorical',xlabels);

    %coulordata=[ones(1,15) 2*ones(1,15) ones(1,14) 2*ones(1,14) ones(1,15) 2*ones(1,15) ones(1,9) 2*ones(1,9)];
    coulordata=[ones(1,15) 2*ones(1,15) ones(1,14) 2*ones(1,14)];

    boxchart(positionaldata, x,"GroupByColor",coulordata)
    grid on
    title(metrics_labels(metric),'interpreter', 'none','FontSize',20,'FontWeight','normal','FontName','Arial')
    xline(1.5)

    legend(["N1", "N2",""])
    set(gca,'FontSize',20)

end

clear metric xlabels coulordata x xgroup i group norm positionaldata







%% Analysis of results
metrics_labels=get_label_metrics(3,node_labels);

comparisons=[1 2];

ttest_results=cell(size(allmetrics_v2));
for i=1:length(allmetrics_v2)
    metrics=allmetrics_v2{i};
    ttest_results{i} = ttest_compare_v2(metrics,metrics_labels,version,length(node_labels),comparisons);
end
%writetable(ttest_results, 'ttest_results.xlsx');
clear comparisons comparison_HCvsP comparison_MRtrixvsFSL comparison_cycle i

%% For visualization in BrainNet nodes AAL116
% nodal metrics:
norm=1;

nodestrength=(349:464); bc=(1:116); lC=(117:232); ec=(233:348);
m=[nodestrength;bc;lC;ec];
names=["D" "BC" "Ci" "EC"];
ttest_results2=ttest_results{norm};
for i=1:1
    qvalues=table2array(ttest_results2(m(i,:),5)); 
    diff=table2array(ttest_results2(m(i,:),7));
    nodes_degree_color = nodes_color_size(qvalues,diff,0.05);
    nodefile = table(makenodefile("aal116_MNIcoord.txt",node_labels,nodes_degree_color));
    writetable(nodefile, 'nodes/ismrm23/'+names(i)+'_midinter_n'+string(norm)+'.txt','Delimiter',' ','WriteVariableNames', 0);
end
clear pvalues diff nodes_degree_color nodefile nodestrength bc lC ec i qvalues pvals m diff norm

%%
names=["BC" "Ci" "EC" "D"];
node_labels=get_label_nodes("AAL116_labels_number.txt");
for i=1:4
    qvalues=table1(:,i); 
    diff=ones(116,1);
    nodes_degree_color = nodes_color_size(qvalues,diff,0.05/116/4);
    nodefile = table(makenodefile("aal116_MNIcoord.txt",node_labels,nodes_degree_color));
    writetable(nodefile, 'nodes/ismrm23/'+names(i)+'_diffnorms.txt','Delimiter',' ','WriteVariableNames', 0);
end



%% Plot histogram ROI sizes

dir_roi='/Users/ana/Documents/Ana/universidade/Tese/Code/matlab_scripts/roi_sizes_intersect';
F = dir_roi;
filePattern = fullfile(F,"*_intersect*");
theFiles = dir(filePattern);
roisizes=zeros(116,1);
for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(theFiles(k).folder, baseFileName);
    %disp(fullFileName)
    roi_size=importdata(fullFileName); 
    if contains(fullFileName,"midcycle")
        idx=1;
    elseif contains(fullFileName,"interictal")
        idx=2;
    elseif contains(fullFileName,"premenstrual")
        idx=3;
        continue
    else
        idx=4;
        continue
    end

    roisizes=horzcat(roisizes,roi_size);
end
%c=importdata("aal.csv");
figure('color','w','Position',[386 100 895 500])
set(gca,'FontSize',20)
bar(mean(roisizes, 2))
xticks([0 10 30 50 70 90 110])
%set(gca,'XTickLabel',c);
%rotateXLabels(gca, 45);
xlabel("Region index in AAL116")
ylabel("Average Region Volume (in voxels)")


clear x roi_size idx fullFileName baseFileName k theFiles filePattern F dir_roi

