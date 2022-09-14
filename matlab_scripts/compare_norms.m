%% Compare Normalizations/Algorithms/conditions
clc
clear variables
close all
format long
%% Load data from matrices

dir='/Users/ana/Documents/Ana/universidade/Tese/Code/matlab_scripts/matrix_data';
dir_roi='/Users/ana/Documents/Ana/universidade/Tese/Code/matlab_scripts/roi_sizes';
atlas="AAL116";
threshold=0;
normalizations=[1 2 3 4];
[allconnectomes,n_conditions,n_people,node_labels,condition_names] = get_data(dir,dir_roi,atlas,threshold,normalizations,false);

clear threshold atlas dir dir_roi
%% Scatter plot and correlation coefficients

nnodes=length(node_labels);
for con=1:length(allconnectomes)
    scatterv=[[],[]];
    scatterv=zeros((sum(n_people)/2)*((nnodes*nnodes-nnodes)/2),2);
    idx=1;
    connectomes=allconnectomes{con};
    for c=[1 2 5 6]%linspace(1,7,4)
        for p=1:n_people(c)
            for i=1:nnodes-1
                for j=i+1:nnodes
                    scatterv(idx,:)=[connectomes{c}(i,j,p) connectomes{c+2}(i,j,p)];
                    idx=idx+1;
                end
            end
        end
    end

    scattervlog=log10(scatterv);
    for i=length(scattervlog):-1:1
        if isinf(scattervlog(i,1)) || isinf(scattervlog(i,2))
            scattervlog(i,:)=[];
        end
    end

    figure('color','w')
    scatterhist(scatterv(:,1),scatterv(:,2), 'NBins',[40,40],'Direction','out',Marker='x');
    xlabel('MRTrix');ylabel('FSL')
    title("Normal Scale: Normalization "+num2str(con))


    figure('color','w')
    scatterhist(scattervlog(:,1),scattervlog(:,2), 'NBins',[40,40],'Direction','out',Marker='x');
    xlabel('MRTrix');ylabel('FSL')
    title("Log Scale: Normalization "+num2str(con))

    R= corr(scatterv,"Type","Pearson");
    disp("Normal Scale Pearson's correlation coefficient N"+num2str(con)+": "+ num2str(R(1,2)))
    Rlog= corr(scattervlog,"Type","Pearson");
    disp("Log Scale Pearson's correlation coefficient N"+num2str(con)+": "+ num2str(Rlog(1,2)))
end

clear i j p c R Rlog scatterv idx con
%% Get metrics


allmetrics=cell(size(allconnectomes));
version_metrics=2;%  3=nodal metrics, 2=general metrics
load('allmetrics2.mat')
metrics_labels=get_label_metrics(version_metrics,node_labels);

for i=1:length(allconnectomes)
    connectomes=allconnectomes{i};
    allmetrics{i}=get_metrics(connectomes,version_metrics);
end

clear connectomes i version_metrics


%% Do kruskall wallis - Compare Groups
n_norms=length(allmetrics);

cond=[1 2 1 2 3 4 3 4];
g_alg=[1 1 2 2 1 1 2 2];
txt = input("Do you want to use multcompare?[y/n]");

for metric=1:length(metrics_labels) % for all metrics
    for norm=1:n_norms % for all normalizations
        metrics_norm=allmetrics{norm};
        for alg=1:2 % for all algorithms
            i=1;
            data=[];
            g_cond=[];
            for conds=1:numel(metrics_norm)% for all conditions
                if g_alg(conds)==alg
                    metrics=cell2mat(metrics_norm(conds));
                    for data_point=1:length(metrics(metric,:))
                        data=[data metrics(metric,data_point)];
                        g_cond=[g_cond cond(conds)];
                        i=i+1;
                    end
                end
            end
            [p_group,~,stats] = kruskalwallis(data,g_cond,"off");
            if p_group <= 0.05
                disp("p-value "+metrics_labels(metric)+" for N"+norm+", A"+alg+ ": "+p_group)
                if txt=='y'
                    figure;
                    multcompare(stats)
                else %txt=='n'
                    for cond1=1:numel(metrics_norm)-1
                        for cond2=cond1+1:numel(metrics_norm)
                            [x,y] = extract_groups(data,g_cond,cond1,cond2);
                            if all([cond1 cond2]==[1 2]) || all([cond1 cond2]==[3 4]) % HC vs M->wilcoxon ranksum
                                p=ranksum(x,y);
                                disp("p-value "+cond1+"-"+cond2+ ": "+p)
                            elseif length(x)==length(y) && (all([cond1 cond2]==[1 3]) || all([cond1 cond2]==[2 4])) % Cycle->wilcoxon signed rank
                                p=signrank(x,y);
                                disp("p-value "+cond1+"-"+cond2+ ": "+p)
                            end
                        end
                    end
                end
            end
        end
    end
end

clear cond data i g_cond metric norm conds metrics_norm metrics data_point p stats tbl x y txt cond1 cond2 g_alg

%% Do Friedman test - Compare Normalization (for all conditions and algorithms)

n_norms=length(allmetrics);

g_alg=[1 1 2 2 1 1 2 2];
reps=sum(n_people)/2;
txt = input("Do you want to use multcompare?[y/n]");
Table =zeros (length(metrics_labels)*sum(n_people),n_norms)

%T = table(LastName,Gender,Age,Height,Weight,Smoker,Systolic,Diastolic);
for metric=1:length(metrics_labels) % for all metrics
    data1=zeros(sum(n_people),n_norms);
    data2=zeros(sum(n_people),n_norms);
    for norm=1:n_norms % for all normalizations
        i=1;j=1;
        metrics_norm=allmetrics{norm};
        for alg=1:2 % for all algorithms
            for conds=1:numel(metrics_norm)% for all conditions
                metrics=cell2mat(metrics_norm(conds));
                for data_point=1:length(metrics(metric,:))
                    if alg==g_alg(conds)
                        data1(i,norm)=metrics(metric,data_point);
                        i=i+1;
                    end
                    if alg==1
                        data2(j,norm)=metrics(metric,data_point);
                    end
                    j=j+1;
                end
            end
        end
    end
    disp(" ")
    [p,~,stats]=friedman(data1,reps,"off");

    if p<0.05
        disp("p-value "+metrics_labels(metric)+ " : "+p)

        if txt=='y'
            multcompare(stats, "Display","off")
        else
            for cond1=1:numel(metrics_norm)
                for norm1=1:3
                    for norm2=norm1+1:4
                        idx=1+sum(n_people(1:cond1))-n_people(cond1):sum(n_people(1:cond1));
                        x=data2(idx,norm1);
                        y=data2(idx,norm2);
                        p_signrank=signrank(x,y);
                        if p_signrank>0.05
                            disp("p-value "+cond1+", "+norm1 +"-"+norm2+": "+p_signrank)
                        end
                    end
                end
            end
        end
    end
end

%% Do Friedman test - Compare Algorithm (for all normalizations and conditions and for each normalization)

n_norms=length(allmetrics);

g_alg=[1 1 2 2 1 1 2 2];
reps=sum(n_people)/2;
txt = input("Do you want to account for all normalizations?[y/n]");
for metric=1:length(metrics_labels) % for all metrics
    data=zeros(sum(n_people)/2*n_norms,2);
    for alg=1:2 % for all normalizations
        i=1;
        for norm=1:n_norms
            metrics_norm=allmetrics{norm};
            for conds=1:numel(metrics_norm)% for all conditions
                if g_alg(conds)==alg
                    metrics=cell2mat(metrics_norm(conds));
                    for data_point=1:length(metrics(metric,:))
                        data(i,alg)=metrics(metric,data_point);
                        i=i+1;
                    end
                end
            end
        end
    end

    
    if txt=='y'
        [p,tbl,stats]=friedman(data,reps,"off");
        if p<0.05
            disp("p-value "+metrics_labels(metric)+ " : "+p)
            multcompare(stats,"Display","off")
        end
    else % Separate normalizations
        for i=1:n_norms
            [p,~,stats]=signrank(data((i-1)*reps+1:reps*i,1),data((i-1)*reps+1:reps*i,2));
            if p<0.05
                disp("p-value "+metrics_labels(metric)+" for N"+i+ " : "+p)
            end
        end
    end 
end
clear txt p norm i conds alg data_point metric reps
%% Do Friedman test - Compare Cycle stages
n_norms=length(allmetrics);

cond=[1 2 1 2 3 4 3 4];
g_alg=[1 1 2 2 1 1 2 2];

for metric=1:length(metrics_labels) % for all metrics
    for norm=1:n_norms % for all normalizations
        metrics_norm=allmetrics{norm};
        for alg=1:2 % for all algorithms
            i=1;
            data=[];
            g_cond=[];
            for conds=1:numel(metrics_norm)% for all conditions
                if g_alg(conds)==alg
                    metrics=cell2mat(metrics_norm(conds));
                    for data_point=1:length(metrics(metric,:))
                        data=[data metrics(metric,data_point)];
                        g_cond=[g_cond cond(conds)];
                        i=i+1;
                    end
                end
            end
            [p,~,stats] = kruskalwallis(data,g_cond,"off");
            if p<0.05
                disp("p-value "+metrics_labels(metric)+" for N"+norm+", A"+alg+ ": "+p)
                figure;
                multcompare(stats)
            end
        end
    end
end

clear cond data i g_cond metric norm conds metrics_norm metrics data_point p stats tbl

%% Do ANOVA
alg=[1 1 2 2 1 1 2 2];
cond=[1 2 1 2 3 4 3 4];

data=zeros(1,n_norms*sum(n_people));
g_alg=zeros(1,n_norms*sum(n_people));
g_norm=zeros(1,n_norms*sum(n_people));
g_cond=zeros(1,n_norms*sum(n_people));
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
    
end

clear i data_point conds norm metric metrics_norm metrics alg cond
%% Plot Boxchart (color=algorithm)

notlog=[3 4 6 7];

groups_color=[1 1 2 2 1 1 2 2];


for metric=1:7
    clear data;i=1;
    for norm=1:length(normalizations)
        metrics=allmetrics{norm};
        for c=1:n_conditions
            data{i,:}=metrics{c}(metric,:);
            i=i+1;
        end
    end
    data=data';
    positionaldata=[];
    colordata=[];
    groups_position=[1 2 1 2 3 4 3 4 5 6 5 6 7 8 7 8 9 10 9 10 11 12 11 12 13 14 13 14 15 16 15 16];
    groups_color=repmat(groups_color,[1 length(normalizations)]);


    for i=1:length(data)
        positionaldata=[positionaldata groups_position(i)*ones(1,length(data{i}))];
        colordata=[colordata groups_color(i)*ones(1,length(data{i}))];

    end

    %     disp(metrics_labels(metric))
    %     for i=1:length(data)-1
    %         for j=i+1:length(data)
    %             p=stattest(data{i},data{j});
    %             if p<0.05/(16*15)
    %                 disp(num2str(i)+"-"+num2str(j)+": p="+num2str(p))
    %             end
    %         end
    %     end

    data=cell2mat(data);
    colorlabels={'MRTrix' 'FSL'};
    xlabels={'HCmid-N1' 'Mint-N1' 'HCpre-N1' 'Mict-N1' ...
        'HCmid-N2' 'Mint-N2' 'HCpre-N2' 'Mict-N2'...
        'HCmid-N3' 'Mint-N3' 'HCpre-N3' 'Mict-N3'...
        'HCmid-N4' 'Mint-N4' 'HCpre-N4' 'Mict-N4'};
    positionaldata=discretize(positionaldata,1:17,'categorical',xlabels);

    figure('color','w','units','normalized','Position',[0.2,0.2,0.8,0.5])
    boxchart(positionaldata,data,'GroupByColor',colordata);
    if ~ismember(metric,notlog)
        set(gca, 'YScale', 'log');
    end
    grid on
    legend(colorlabels)
    title(metrics_labels(metric),'interpreter', 'none')
end

clear notlog positionaldata xlabels colorlabels data colordata groups_color groups_position i c metric norm

%% Plot Boxchart (color=condition)

notlog=[3 4 6 7];

for metric=1:7
    clear data;i=1;
    for norm=1:length(normalizations)
        metrics=allmetrics{norm};
        for c=1:n_conditions
            data{i,:}=metrics{c}(metric,:);
            i=i+1;
        end
    end
    data=data';
    positionaldata=[];
    colordata=[];
    groups_position=[1 1 2 2 1 1 2 2 3 3 4 4 3 3 4 4 5 5 6 6 5 5 6 6 7 7 8 8 7 7 8 8]; %algorithm and norm
    groups_color=[1 2 1 2 3 4 3 4 1 2 1 2 3 4 3 4 1 2 1 2 3 4 3 4 1 2 1 2 3 4 3 4]; %conditions

    for i=1:length(data)
        positionaldata=[positionaldata groups_position(i)*ones(1,length(data{i}))];
        colordata=[colordata groups_color(i)*ones(1,length(data{i}))];

    end

    data=cell2mat(data);
    colorlabels={'HCmid' 'Minter' 'HCpre' 'Mict'};
    xlabels={'MRtrix-N1' 'FSL-N1' ...
        'MRtrix-N2' 'FSL-N2' ...
        'MRtrix-N3' 'FSL-N3'...
        'MRtrix-N4' 'FSL-N4'};
    positionaldata=discretize(positionaldata,1:9,'categorical',xlabels);

    figure('color','w','units','normalized','Position',[0.2,0.2,0.8,0.5])
    boxchart(positionaldata,data,'GroupByColor',colordata);
    if ~ismember(metric,notlog)
        set(gca, 'YScale', 'log');
    end
    grid on
    legend(colorlabels)
    title(metrics_labels(metric),'interpreter', 'none')
end

clear notlog positionaldata xlabels colorlabels data colordata groups_color groups_position i c metric norm