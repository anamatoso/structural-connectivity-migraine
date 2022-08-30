%% Compare Normalizations

%% Load data from matrices
clear variables
close all
format long
dir='/Users/ana/Documents/Ana/universidade/Tese/Code/matlab_scripts/matrix_data';
dir_roi='/Users/ana/Documents/Ana/universidade/Tese/Code/matlab_scripts/roi_sizes';
atlas="AAL116"; 
threshold=0;
normalizations=[1 3];
[allconnectomes,n_conditions,n_people,node_labels,condition_names] = get_data(dir,dir_roi,atlas,threshold,normalizations,false);

clear threshold atlas dir dir_roi
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





