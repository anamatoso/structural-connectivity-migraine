function [connectomes] = load_data(folder_matrix,pattern,folder_roi_sizes,algorithm)
% This function loads the data from a given folder which have a given
% pattern

F = folder_matrix;
filePattern = fullfile(F, pattern);
theFiles = dir(filePattern);
if isempty(theFiles)
    disp ('Error: File of matrices is empty')
    return 
end
baseFileName = theFiles(1).name;
fullFileName = fullfile(theFiles(1).folder, baseFileName);
c=importdata(fullFileName);
connectomes=zeros(length(c),length(c),length(theFiles));

for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    if algorithm=="fsl"
        roi_file=strcat(folder_roi_sizes,"/",baseFileName(1:length(baseFileName)-17),"_roi_size.txt"); %17
        fullFileName = fullfile(theFiles(k).folder, baseFileName);
        c=normalize_roisize_fsl(fullFileName,roi_file);
        
    elseif algorithm=="mrtrix"
        roi_file=strcat(folder_roi_sizes,"/",baseFileName(1:length(baseFileName)-34),"_roi_size_intersect.txt"); %17
        fullFileName = fullfile(theFiles(k).folder, baseFileName);
        c=normalize_roisize_mrtrix(fullFileName,roi_file);
    
    else
        disp ('Error: Algorithm not known. Please choose either fsl or mrtrix.')
        return
    end
    connectomes(:,:,k)=c;
end
end
