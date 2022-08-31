%% Compare Normalizations
clear variables
close all
format long
%% Load data from matrices

dir='/Users/ana/Documents/Ana/universidade/Tese/Code/matlab_scripts/matrix_data';
dir_roi='/Users/ana/Documents/Ana/universidade/Tese/Code/matlab_scripts/roi_sizes';
atlas="AAL116"; 
threshold=0;
normalizations=[1 2];
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
    scatterhist(scattervlog(:,1),scattervlog(:,2), 'NBins',[40,40],'Direction','out',Marker='x'); 
    xlabel('MRTrix');ylabel('FSL')
    title("Normalization "+num2str(con))
    R= corr(scattervlog,"Type","Pearson");
    disp(" ")
    disp("Pearson's correlation coefficient N"+num2str(con)+": "+ num2str(R(1,2)))
end


%% Get metrics
allmetrics=cell(size(allconnectomes));
version_metrics=2;%  1=nodal metrics, 2=general metrics

for i=1:length(allconnectomes)
    connectomes=allconnectomes{i};
    allmetrics{i}=get_metrics(connectomes,version_metrics);
end
metrics_labels=get_label_metrics(version_metrics,node_labels);

clear connectomes i version_metrics

%% Plot Boxchart

notlog=[3 4 6 7];
for metric=1:7
    clear data;i=1;
    for norm=1:length(normalizations)
        metrics=allmetrics{normalizations(norm)};
        for c=1:n_conditions
            data{i,:}=metrics{c}(metric,:);
            i=i+1;
        end
    end
    data=data';
    positionaldata=[];
    colordata=[];
    groups_position=[1 2 1 2 3 4 3 4 5 6 5 6 7 8 7 8];
    groups_color=[1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2];
    
    for i=1:length(data)
        positionaldata=[positionaldata groups_position(i)*ones(1,length(data{i}))];
        colordata=[colordata groups_color(i)*ones(1,length(data{i}))];
        
    end

    disp(metrics_labels(metric))
    for i=1:length(data)-1
        for j=i+1:length(data)
            p=stattest(data{i},data{j});
            if p<0.05/(16*15)
                disp(num2str(i)+"-"+num2str(j)+": p="+num2str(p))
            end
        end
    end



    data=cell2mat(data);
    colorlabels={'MRTrix' 'FSL'};
    xlabels={'HCmid-N1' 'Mint-N1' 'HCpre-N1' 'Mict-N1' 'HCmid-N2' 'Mint-N2' 'HCpre-N2' 'Mict-N2'};
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





