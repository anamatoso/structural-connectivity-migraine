function [allconnectomes,n_conditions,n_people,node_labels,condition_names] = get_data(dir,dir_roi,atlas,threshold,normalizations)
% This function gets the data of all connectivity matrices and saves them in allconnectomes variable. It also gets the number of people and the number of groups. Additionally it also creates the node labels and the labels of the groups.

% Decide what is the pattern to look for in the files
if atlas=="AAL116" 
    pattern="_prob";%"_intersect"; 
else
    pattern="*"+atlas; 
end

% Create the allconnectomes variable
i=1;
allconnectomes=cell(length(normalizations),1);
for norm=normalizations
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

    allconnectomes{i}={HC_midcycle_mrtrix HC_midcycle_fsl HC_premenstrual_mrtrix HC_premenstrual_fsl;...
    M_interictal_mrtrix M_interictal_fsl M_ictal_mrtrix M_ictal_fsl};
    i=i+1;
end

% Get number of groups
connectomes=allconnectomes{1};
n_conditions=numel(connectomes);

% Calculate people per group
n_people=zeros(1,n_conditions);
for i=1:n_conditions
    conmat=connectomes{i};
    s=size(conmat);
    n_people(i)=s(end);
end

% Create node labels and the labels of the groups.
node_labels=get_label_nodes(atlas+"_labels.txt");
condition_names=["MRtrix-HC-midcycle" "MRtrix-M-interictal" "FSL-HC-midcycle" "FSL-M-interictal" "MRtrix-HC-premenstrual" "MRtrix-M-ictal" "FSL-HC-premenstrual" "FSL-M-ictal"];

end