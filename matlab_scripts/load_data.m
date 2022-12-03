function [connectomes] = load_data(folder_matrix,pattern,folder_roi_sizes,algorithm,threshold,normalization)
% This function loads the data from a given folder which have a given
% pattern

% Get folder and files
F = folder_matrix;
filePattern = fullfile(F, pattern);
theFiles = dir(filePattern);

% Check if folder is empty
if isempty(theFiles)
    error('Error: File of matrices is empty') 
end

% Import data
baseFileName = theFiles(1).name;
fullFileName = fullfile(theFiles(1).folder, baseFileName);
c=load(fullFileName);
connectomes=zeros(length(c),length(c),length(theFiles)); % create matrix of size [n_nodes,n_nodes,n_people_of_group]

for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(theFiles(k).folder, baseFileName);

    if algorithm=="fsl"
        roi_file=strcat(folder_roi_sizes,"/",baseFileName(1:length(baseFileName)-23),"_roi_size_intersect.txt"); %17, or 23 for omat3
    elseif algorithm=="mrtrix"
        roi_file=strcat(folder_roi_sizes,"/",baseFileName(1:length(baseFileName)-29),"_roi_size_intersect.txt"); %34 aal, 32 desikan
        %disp(roi_file)
    else
        error('Error: Algorithm not known. Please choose either fsl or mrtrix.')
        return
    end
    c=normalize_matrix(fullFileName,roi_file,threshold,algorithm,normalization); %normalize matrix according to certain normalization
    connectomes(:,:,k)=c;
end
end

