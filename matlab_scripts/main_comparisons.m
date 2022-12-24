%% Compare Normalizations/Algorithms/Groups
clc
clear variables
close all
format long

%% Load data from matrices
dir=strcat(pwd,'/matrix_data_prob');
dir_roi=strcat(pwd,'/roi_sizes');
atlas="AAL116";
threshold=0; % threshold of minimum connections
normalizations=[1 2 3 4];
[allconnectomes,n_conditions,n_people,node_labels,condition_names] = get_data(dir,dir_roi,atlas,threshold,normalizations);

clear threshold atlas dir dir_roi

%% Scatter plot and correlation coefficients

n_nodes=length(node_labels);
for con=1:length(allconnectomes)
    scatterv=zeros((sum(n_people)/2)*((n_nodes*n_nodes-n_nodes)/2),2);
    idx=1;
    connectomes=allconnectomes{con};
    for c=[1 2 5 6] %MRtrix connectomes
        for p=1:n_people(c)
            for i=1:n_nodes-1
                for j=i+1:n_nodes
                    scatterv(idx,:)=[connectomes{c}(i,j,p) connectomes{c+2}(i,j,p)];
                    idx=idx+1;
                end
            end
        end
    end

    % Logarithmic relationship
    scattervlog=log10(scatterv);
    for i=length(scattervlog):-1:1
        if isinf(scattervlog(i,1)) || isinf(scattervlog(i,2))
            scattervlog(i,:)=[];
        end
    end

    % Plot linear scale
    figure('color','w')
    scatterhist(scatterv(:,1),scatterv(:,2), 'NBins',[40,40],'Direction','out',Marker='x');
    xlabel('MRTrix');ylabel('FSL')
    title("Linear Scale: Normalization "+num2str(con))

    % Plot logarithmic scale
    figure('color','w')
    scatterhist(scattervlog(:,1),scattervlog(:,2), 'NBins',[40,40],'Direction','out',Marker='x');
    xlabel('MRTrix');ylabel('FSL')
    title("Log Scale: Normalization "+num2str(con))

    % Calculate correlation coefficient
    R= corr(scatterv,"Type","Pearson");
    disp("Normal Scale Pearson's correlation coefficient N"+num2str(con)+": "+ num2str(R(1,2)))
    Rlog= corr(scattervlog,"Type","Pearson");
    disp("Log Scale Pearson's correlation coefficient N"+num2str(con)+": "+ num2str(Rlog(1,2)))
end

clear i j p c R Rlog scatterv idx con scattervlog

%% Plot histogram ROI sizes

% Get Roi sizes' files and sort them according to their group
dir_roi=strcat(pwd,'/roi_sizes_intersect');
F = dir_roi;
filePattern = fullfile(F,"*_intersect*");
theFiles = dir(filePattern);
roisizes=cell(1,4);
for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(theFiles(k).folder, baseFileName);
    
    roi_size=load(fullFileName); 
    if contains(fullFileName,"midcycle")
        idx=1;
    elseif contains(fullFileName,"interictal")
        idx=2;
    elseif contains(fullFileName,"premenstrual")
        idx=3;
    else
        idx=4;
    end

    roisizes{idx}=vertcat(roisizes{idx},roi_size);
end

%titles
condition={'Controls in the Midcycle Phase', 'Migraineurs in the Interictal Phase', ...
    'Controls in the Prementrual Phase', 'Migraineurs in the Ictal Phase'};

% Plot the 4 histograms (1 for each group)
for k=1:4
    figure('color','w','Position',[268 165 652 532])
    histogram(roisizes{k},'Normalization','probability')
    ylabel('Probability')
    title("ROI sizes of "+ condition(k),'interpreter', 'none','FontSize',15,'FontWeight','normal','FontName','Arial')
    set(gca,'FontSize',15)
end

clear roi_size idx fullFileName baseFileName k theFiles filePattern F dir_roi

%% Calculate metrics

allmetrics=cell(size(allconnectomes));
version_metrics=2;%  3=nodal metrics, 2=general metrics
metrics_labels=get_label_metrics(version_metrics,node_labels);

for i=1:length(allconnectomes)
    connectomes=allconnectomes{i};
    allmetrics{i}=get_metrics(connectomes,version_metrics);
end

clear connectomes i version_metrics

%% Do kruskall wallis - Compare Groups (match)

n_norms=length(allmetrics);
total_people=sum(n_people);

cond=[1 2 1 2 3 4 3 4];
g_alg=[1 1 2 2 1 1 2 2];

txt = input("Do you want to use multcompare?['y'/'n']");

for metric=1:length(metrics_labels) % for all metrics
    for norm=1:n_norms % for all normalizations

        % Get data
        metrics_norm=allmetrics{norm};
        for alg=1:2 % for all algorithms
            i=1;
            data=zeros(1,total_people/2);
            g_cond=zeros(1,total_people/2);
            for conds=1:numel(metrics_norm)% for all conditions
                if g_alg(conds)==alg
                    metrics=cell2mat(metrics_norm(conds));
                    for data_point=1:length(metrics(metric,:))
                        data(i)=metrics(metric,data_point);
                        g_cond(i)=cond(conds);
                        i=i+1;
                    end
                end
            end

            % Perform Kruskal-Wallis test
            [p_group,~,stats] = kruskalwallis(data,g_cond,"off");

            % Perform post-hoc tests
            if p_group <= 0.05
                disp("p-value "+metrics_labels(metric)+" for N"+norm+", A"+alg+ ": "+p_group)
                if txt=='y'
                    multcompare(stats,"Display","off")
                else %txt=='n'
                    for cond1=1:numel(metrics_norm)-1
                        for cond2=cond1+1:numel(metrics_norm)
                            [x,y] = extract_groups(data,g_cond,cond1,cond2);
                            if all([cond1 cond2]==[1 2]) || all([cond1 cond2]==[3 4]) % HC vs M->wilcoxon ranksum
                                p=ranksum(x,y);
                                if p<=0.05
                                    disp("p-value "+cond1+"-"+cond2+ ": "+p)
                                end
                            elseif length(x)==length(y) && (all([cond1 cond2]==[1 3]) || all([cond1 cond2]==[2 4])) % Cycle->wilcoxon signed rank
                                p=signrank(x,y);
                                if p<=0.05
                                    disp("p-value "+cond1+"-"+cond2+ ": "+p)
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

clear cond i alg p_group data metric g_cond norm conds metrics_norm metrics data_point p stats x y txt cond1 cond2 g_alg total_people

%% Do kruskall wallis - Compare Groups (not matched)

n_norms=length(allmetrics);
total_people=sum(n_people);

cond=[1 2 1 2 3 4 3 4];
g_alg=[1 1 2 2 1 1 2 2];

txt = input("Do you want to use multcompare?['y'/'n']");

for metric=1:length(metrics_labels) % for all metrics

    % Get data
    for norm=1:n_norms % for all normalizations
        metrics_norm=allmetrics{norm};
        for alg=1:2 % for all algorithms
            i=1;
            data=zeros(1,total_people/2);
            g_cond=zeros(1,total_people/2);
            for conds=1:numel(metrics_norm)% for all conditions
                if g_alg(conds)==alg
                    metrics=cell2mat(metrics_norm(conds));
                    for data_point=1:length(metrics(metric,:))
                        data(i)=metrics(metric,data_point);
                        g_cond(i)=cond(conds);
                        i=i+1;
                    end
                end
            end

            % Perform Kruskal-Wallis test
            [p_group,~,stats] = kruskalwallis(data,g_cond,"off");

            % Perform post-hoc tests
            if p_group <= 0.05
                disp("p-value "+metrics_labels(metric)+" for N"+norm+", A"+alg+ ": "+p_group)
                if txt=='y'
                    multcompare(stats,"Display","off")
                else %txt=='n'
                    for cond1=1:4-1
                        for cond2=cond1+1:4
                            [x,y] = extract_groups(data,g_cond,cond1,cond2);
                            p=ranksum(x,y);
                            if p<=0.05/4
                                disp("p-value "+cond1+"-"+cond2+ ": "+p)
                            end
                        end
                    end
                end
            end
        end
    end
end

clear cond i alg data g_cond p_group metric norm conds metrics_norm metrics data_point p stats tbl x y txt cond1 cond2 g_alg total_people

%% Do Friedman test - Compare Normalization (for all conditions and algorithms)

n_norms=length(allmetrics);

g_alg=[1 1 2 2 1 1 2 2];
reps=sum(n_people)/2;
txt = input("Do you want to use multcompare?[y/n]");
Table=zeros(length(metrics_labels)*length(n_people)*6,6);
t=1;
data_all=cell(1,length(metrics_labels));
for metric=1:length(metrics_labels) % for all metrics

    % Get data
    data_friedman=zeros(sum(n_people),n_norms);
    data_table=zeros(sum(n_people),n_norms);
    for norm=1:n_norms % for all normalizations
        i=1;j=1;
        metrics_norm=allmetrics{norm};
        for alg=1:2 % for all algorithms
            for conds=1:n_conditions% for all conditions
                metrics=cell2mat(metrics_norm(conds));
                for data_point=1:length(metrics(metric,:))
                    if alg==g_alg(conds)
                        data_friedman(i,norm)=metrics(metric,data_point);
                        i=i+1;
                    end
                    if alg==1
                        data_table(j,norm)=metrics(metric,data_point);
                        j=j+1;
                    end    
                end

            end
        end
    end
    
    % Perform Friedman test
    disp(" ")
    [p,~,stats]=friedman(data_friedman,reps,"off");
    data_all{metric}=data_friedman;
    disp("p-value "+metrics_labels(metric)+ " : "+p)

    % Perform Post-hot tests
    if p<0.05
        if txt=='y'
            multcompare(stats, "Display","off")
        else
            for norm1=1:3
                    for norm2=norm1+1:4
                        p=signrank(data_friedman(:,norm1),data_friedman(:,norm2));
                        %if p<0.05/6/length(metrics_labels)
                            disp("posthoc test,"+ "N"+norm1 +"-N"+norm2+": "+p)
                        %end
                    end
            end
           
            for cond1=1:n_conditions
                for norm1=1:3
                    for norm2=norm1+1:4
                        idx=1+sum(n_people(1:cond1))-n_people(cond1):sum(n_people(1:cond1));
                        x=data_table(idx,norm1);
                        y=data_table(idx,norm2);
                        p_signrank=signrank(x,y);
                        if p_signrank<0.05/6/8/length(metrics_labels) % correct for comparisons of norms and conditions
                            %disp("p-value cond:"+cond1+", N"+norm1 +"-N"+norm2+": "+p_signrank)
                        end
                    end
                end
            end
        end
    end
    for cond1=1:n_conditions
        for norm1=1:3
            for norm2=norm1+1:4
                idx=1+sum(n_people(1:cond1))-n_people(cond1):sum(n_people(1:cond1));
                x=data_table(idx,norm1);
                y=data_table(idx,norm2);
                p_signrank=signrank(x,y);
                Table(t,:)=[metric cond1 norm1 norm2 p_signrank median(y)>median(x)];t=t+1;
            end
        end
    end
end

names_table=cell(length(Table),2);
for i=1:length(Table)
    names_table{i,1}=metrics_labels(Table(i,1));
    names_table{i,2}=condition_names(Table(i,2));
end

table=array2table(Table(:,3:end),'VariableNames',{'Norm1','Norm2','pvalue', 'Variation'});
names_table=cell2table(names_table,'VariableNames',{'Metric','Condition'});
final=[names_table table];


clear alg metrics metrics_norm names_table cond1 conds data_point g_alg i j idx t metric n_norms norm reps x y txt stats p p_signrank norm1 norm2 table Table data_friedman data_table data

%% Do Friedman test - Compare Algorithm (for all normalizations and conditions and for each normalization)

n_norms=length(allmetrics);
n_people2=[15 14 15 9];
reps=sum(n_people)/2;
conds=["Midcycle" "Interictal" "Premenstrual" "Ictal"];
txt = input("Do you want to use multcompare?[y/n]");
Table=zeros(length(metrics_labels)*length(n_people),5);
t=1;

for metric=1:length(metrics_labels) % for all metrics
    
    % Get data
    data_friedman=data_all{metric};    
    disp(" ")

    data1=zeros(reps*n_norms,2);
    for norm=1:n_norms
        data1((norm-1)*reps+1:norm*reps,:)=reshape(data_friedman(:,norm),[reps,2]);
    end
    
    % Perform Friedman test
    [p,~,stats]=friedman(data1,reps,"off");
    disp("p-value "+metrics_labels(metric)+ " : "+p)

    % Perform Post-hoc tests
    if txt=='y'
        multcompare(stats, "Display","off")
    else
        p=signrank(data1(:,1),data1(:,2));
        if p<2%0.05/length(metrics_labels)
            disp(metrics_labels(metric) +" signranktest: "+p)
        end
        for norm1=1:4
            for cond1=1:4
                    idx=norm1*(1+sum(n_people2(1:cond1))-n_people2(cond1)):norm1*(sum(n_people2(1:cond1)));
                    x=data1(idx,1);
                    y=data1(idx,2);
                    p_signrank=signrank(x,y);
                    if p_signrank<0.05/6/4/length(metrics_labels) % correct for norms, conds and n_metrics
                        %disp("p-value N"+norm1+", "+conds(cond1)+": "+p_signrank)
                    end                        
            end
        end  
    end

    % Organize table
    for norm1=1:4
        for cond1=1:4
            idx=norm1*(1+sum(n_people2(1:cond1))-n_people2(cond1)):norm1*(sum(n_people2(1:cond1)));
            x=data1(idx,1);
            y=data1(idx,2);
            p_signrank=signrank(x,y);
            Table(t,:)=[metric cond1 norm1 p_signrank median(y)>median(x)];t=t+1;
        end
    end
end

% Create table
names_table=cell(length(Table),2);
for i=1:length(Table)
    names_table{i,1}=metrics_labels(Table(i,1));
    names_table{i,2}=conds(Table(i,2));
end

table=array2table(Table(:,3:end),'VariableNames',{'Normalisation','Pvalue','Variation'});
names_table=cell2table(names_table,'VariableNames',{'Metric','Condition'});
final_table=[names_table table];

clear n_people2 data1 cond1 conds i idx t metric n_norms norm reps x y txt...
    stats p p_signrank norm1 table Table data_table data_friedman names_table

%% Compare algs - Just for N1
n_norms=length(allmetrics);
n_nodes=length(node_labels);

metrics=allmetrics{1};
mrtrix=[1 2 5 6];
pvals=zeros(1,n_nodes*n_norms);
for metric=1:length(metrics_labels) % for all metrics
    data={[];[]};
    for cond=1:n_conditions
        if any(ismember(mrtrix,cond))
            data{1}=horzcat(data{1},metrics{cond}(metric,:));
        else
            data{2}=horzcat(data{2},metrics{cond}(metric,:));
        end
    end
    data=cell2mat(data);
    p=signtest(data(1,:),data(2,:));
    %if p<=0.5/116/4
        disp(metrics_labels(metric) +" : "+p)
    %end
    pvals(metric)=p;

end

% Nodal metrics
nodestrength=(349:464); bc=(1:116); lC=(117:232); ec=(233:348);
m=[nodestrength;bc;lC;ec];
names=["nodestrength" "bc" "lC" "ec"];
diff=ones(1,n_nodes);
node_labels=1:n_nodes;
for i=1:n_norms
    qvalues=pvals(m(i,:)); 
    nodes_degree_color = nodes_color_size(qvalues,diff,0.05/n_nodes/n_norms);
    nodefile = table(makenodefile("aal116_MNIcoord.txt",node_labels,nodes_degree_color));
    writetable(nodefile, 'nodes/model/prob/'+names(i)+'_significantpval.txt','Delimiter',' ','WriteVariableNames', 0);
end

clear metrics cond n_norms n_nodes metric mrtrix p pvals data nodestrength bc lC ec m names diff node_labels i qvalues nodes_degree_color nodefile
%% Perform ANOVA
alg=[1 1 2 2 1 1 2 2];
cond=[1 2 1 2 3 4 3 4];
n_norms=length(allmetrics);
total_people=sum(n_people);

data=zeros(1,n_norms*total_people);
g_alg=zeros(1,n_norms*total_people);
g_norm=zeros(1,n_norms*total_people);
g_cond=zeros(1,n_norms*total_people);

i=1;
for metric=1:7 % for all metrics
    for norm=1:n_norms % for all normalizations
        metrics_norm=allmetrics{norm};
        for conds=1:numel(metrics_norm)% for all conditions
            metrics=cell2mat(metrics_norm(conds));
            for data_point=1:length(metrics(metric,:))
                data(i)=metrics(metric,data_point);
                g_alg(i)=alg(conds);
                g_norm(i)=norm;
                g_cond(i)=cond(conds);
                i=i+1;
            end
        end
    end

    disp(metrics_labels(metric))
    [p,tbl,stats,terms] = anovan(data,{g_alg,g_norm,g_cond},'model', 'full','varnames',{'Alg','Norm','Cond'});
    % Not finished...
end

clear i total_people data_point conds norm metric metrics_norm metrics alg cond data g_alg g_norm g_cond
%% Plot Boxchart (colour=algorithm) - Only Global Metrics!!

notlog=[3 4 6 7];
groups_color=[1 1 2 2 1 1 2 2];

n_norms=length(normalizations);
total_people=sum(n_people);


for metric=1:length(metrics_labels)
    % Get data
    data=cell(1,n_norms*n_conditions);
    i=1;
    for norm=1:length(normalizations)
        metrics=allmetrics{norm};
        for c=1:n_conditions
            data{i}=metrics{c}(metric,:);
            i=i+1;
        end
    end
    
    % Organise grouping
    groups_position=[1 2 1 2 3 4 3 4 5 6 5 6 7 8 7 8 9 10 9 10 11 12 11 12 13 14 13 14 15 16 15 16];
    groups_color=repmat(groups_color,[1 length(normalizations)]);

    positionaldata=zeros(1,total_people*n_norms);
    colordata=zeros(1,total_people*n_norms);
    past=0;
    for i=1:length(data)
        positionaldata(past+1:past+length(data{i}))=groups_position(i)*ones(1,length(data{i}));
        colordata(past+1:past+length(data{i}))=groups_color(i)*ones(1,length(data{i}));
        past=past+length(data{i});
    end

    data=cell2mat(data);

    % Set labels
    colorlabels={'MRTrix' 'FSL'};
    xlabels={'HCmid-N1' 'Mint-N1' 'HCpre-N1' 'Mict-N1' ...
        'HCmid-N2' 'Mint-N2' 'HCpre-N2' 'Mict-N2'...
        'HCmid-N3' 'Mint-N3' 'HCpre-N3' 'Mict-N3'...
        'HCmid-N4' 'Mint-N4' 'HCpre-N4' 'Mict-N4'};
    positionaldata=discretize(positionaldata,1:17,'categorical',xlabels);
    
    % Plot
    figure('color','w','units','normalized','Position',[0.2,0.2,0.8,0.5])
    boxchart(positionaldata,data,'GroupByColor',colordata);
    if ~ismember(metric,notlog)
        set(gca, 'YScale', 'log');
    end
    grid on
    legend(colorlabels)
    title(metrics_labels(metric),'interpreter', 'none')
end

clear notlog past n_norms total_people xlabels positionaldata colorlabels data colordata groups_color groups_position i c metric norm

%% Plot Boxchart (color=condition)

notlog=[3 4 6 7];

n_norms=length(normalizations);
total_people=sum(n_people);

for metric=1:7

    % Get data
    data=cell(1,n_norms*n_conditions);
    i=1;
    for norm=1:n_norms
        metrics=allmetrics{norm};
        for c=1:n_conditions
            data{i}=metrics{c}(metric,:);
            i=i+1;
        end
    end
    
    % Organise grouping
    groups_position=[1 1 2 2 1 1 2 2 3 3 4 4 3 3 4 4 5 5 6 6 5 5 6 6 7 7 8 8 7 7 8 8]; %algorithm and norm
    groups_color=[1 2 1 2 3 4 3 4 1 2 1 2 3 4 3 4 1 2 1 2 3 4 3 4 1 2 1 2 3 4 3 4]; %conditions

    positionaldata=zeros(1,total_people*n_norms);
    colordata=zeros(1,total_people*n_norms);
    past=0;
    for i=1:length(data)
        positionaldata(past+1:past+length(data{i}))=groups_position(i)*ones(1,length(data{i}));
        colordata(past+1:past+length(data{i}))=groups_color(i)*ones(1,length(data{i}));
        past=past+length(data{i});
    end

    data=cell2mat(data);

    % Set labels
    colorlabels={'HCmid' 'Minter' 'HCpre' 'Mict'};
    xlabels={'MRtrix-N1' 'FSL-N1' ...
        'MRtrix-N2' 'FSL-N2' ...
        'MRtrix-N3' 'FSL-N3'...
        'MRtrix-N4' 'FSL-N4'};
    positionaldata=discretize(positionaldata,1:9,'categorical',xlabels);

    % Plot
    figure('color','w','units','normalized','Position',[0.2,0.2,0.8,0.5])
    boxchart(positionaldata,data,'GroupByColor',colordata);
    
    if ~ismember(metric,notlog)
        set(gca, 'YScale', 'log');
    end
    grid on
    legend(colorlabels)
    title(metrics_labels(metric),'interpreter', 'none', 'FontWeight','normal','FontName','Arial','FontSize',20)
end

clear notlog past positionaldata n_norms xlabels colorlabels data colordata groups_color groups_position i c metric norm