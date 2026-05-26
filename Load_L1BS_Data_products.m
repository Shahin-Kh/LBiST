function L1BS=Load_L1BS_Data_products(folderPath)
filePattern = fullfile(folderPath, '**', '*.nc'); 
files = dir(filePattern);
filePaths = {};
L1BS = cell(1, length(filePaths));
    for k = 1:length(files)
        fullPath = fullfile(files(k).folder, files(k).name);
        L1BS{k} = (fullPath);
    end
end
