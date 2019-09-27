function [lfpStruct] = openephys2struct(path, dataDir, saveStruct, metaRat)
% Convert an open eophys set of files to the lfp struct
% requires load_open_ephys_data_faster 
% 
% init 190910 kwc ud 190926 kwc


fprintf('Importing the data from OpenEphys format...\n')
%fNamePath = 'C:\Users\FortinLab\Documents\Data\SG15_sleep\SG15_sleep_S1_2019-09-03_17-33-27\';
%fNamePath = 'C:\Users\FortinLab\Documents\Data\SM04_task\SM04_Session262a_2019-08-14_10-42-17\';

fNamePath = [path, dataDir];
files = dir([fNamePath,'*_CH*.continuous']);

num_channels = length(files);
dsRate = 500;
dsBool = true;
% saveStruct = true;

fprintf('Loading(%i): ',length(files)); tic;

ch = 1; % run on first channel only to get size for allocation
fName = [fNamePath, '100_CH', num2str(ch), '.continuous']; fprintf('%i ',ch);
[tmpData, timestamps, info] = load_open_ephys_data_faster(fName);

% downsample the data
% reduced size from 20GB to .016GB
if dsBool
    fs = info.header.sampleRate;
    dsStep = ceil(fs / dsRate);
    tmpData = tmpData(1:dsStep:end); 
end

% Allocate the matrix for speed
lfpData = nan(num_channels, length(tmpData));
lfpData(ch,:) = tmpData';

for ch = 2:num_channels
    fName = [fNamePath, '100_CH', num2str(ch), '.continuous']; fprintf('%i ',ch);
    [tmpData, timestamps, info] = load_open_ephys_data_faster(fName);
    
    % downsample the data
    if dsBool
    fs = info.header.sampleRate;
    dsStep = ceil(fs / dsRate);
    tmpData = tmpData(1:dsStep:end); 
    end
    
    % store the data
    lfpData(ch,:) = tmpData';
end

if dsBool; fs = fs/dsStep; end

% Wrap it all up with a bow and ribbon
lfpStruct = struct;
lfpStruct.data = lfpData;
lfpStruct.info = info;
lfpStruct.info.dataDir = dataDir; % TODO  change this to the metarat one
lfpStruct.info.dsBool = dsBool;
lfpStruct.info.fs = fs;
lfpStruct.info.date = date;
lfpStruct.info.origFormat = 'OpenEPhys';

lfpStruct.info.ratName = metaRat.ratName; 

% Save our hard work!
if saveStruct; save(lfpStruct.info.dataDir(1:end-1), 'lfpStruct'); 
    fprintf(['\nLFP Saved as ', lfpStruct.info.dataDir(1:end-1), '!\n']); end

fprintf('\nFinished! \n'); toc;

end