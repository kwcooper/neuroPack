function [LFP_Data] = importTrodesLFP(path, dataDir)
% A wrapper script to load in the Trodes data
% Feed it a path to the extracted LFP files; can use a script or you can
% used the trodes software to build this
%

% TODO: 
%   Add metadata
%   fix the .timestamps issue
%   Add verbose optional argument
%   preallocate the data struct; not slow, but could be in the future...

cd(path) % Move to the proper top level dir

% Gather all of the data files in the directory
files = dir([dataDir '*.dat']);

% TODO:
% Be weary of the timestamps file, a better fix will be needed
% but for now just boot the last file, hoping that it will always
% be in the same place...
files(33) = [];

% iterate through the data files and import the data into matlab
% I'm impressed by how quick this is...
fprintf('Reading in the files')
for i = 1:length(files)
    fName = [dataDir, files(i).name];
    % Now extract the data and append it to our larger struct
    data(i) = readTrodesExtractedDataFile(fName);
    
    if mod(i,10) == 0; fprintf('.'); end % Loading bar
end

fprintf('\nfin\n')
LFP_Data = data;

% timestamps file testing
if 0 
    fName = [dataDir, 'SM03_task1_191008.timestamps.dat']; 
    ts_data = readTrodesExtractedDataFile(fName);
end

end


