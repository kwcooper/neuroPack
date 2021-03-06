
function [lfpStruct] = trodes2struct(path, dataDir, saveStruct, metaRat)
% Script to get a look at the trodes data. Imports, downsamples, sorts, and
% saves to a more useful struct for further analysis. 
% Requires importTrodesLFP
% init 190910 kwc ud 191016 kwc

% Call the wrapper function that grabs all of the channels and stores them
% in a handy struct for further use.
% path = 'C:\Users\FortinLab_OE1\Desktop\trodes';
% dataDir = 'SM03-cageTest-190905.LFP/';

%path = 'C:\Users\FortinLab\Documents\Data\SMO3_piloting';
%dataDir = 'SM03-cageTest-190905.LFP/';
fprintf('Importing the data from Trodes format...\n')
LFP_Data = importTrodesLFP(path, dataDir);

%fs = LFP_Data(1).clockrate; % TODO: NEED TO CONFIRM IF THIS IS THE FS!!! (it isn't.. see below)
fs = 1500;        % as per the downsampling from exportlfp and mattia's recomendation
fs_new = 1000;    % desired sampling rate (if dsBool is true)
dsBool = true;    % do ya wanna downsample? (enabled by default)
% saveStruct = true;       
showBuildMat = 0;        % If we want to see how the lfp is being sorted...

% So the trodes import function doesn't keep the data in order when it
% imports it... We can reorder it now, and clean up the struct a lil bit
% iterate through the lfp data struct from trodes and grab and reformat 

% Do some precalculations to figure out sizes
tmpData = LFP_Data(1).fields.data';
    
% downsample the data
if dsBool; tmpData = resampSig(tmpData, fs, fs_new); end

% allocate an empty array for sorting
lfpMat = nan(length(LFP_Data), size(tmpData,2)); 

if showBuildMat; subInd = 1; figure; end

fprintf('Loading(%i): ',length(LFP_Data)); tic;
for ch = 1:length(LFP_Data)
    % This is important! to see whcih channel the lfp belongs to. 
    ntrode = LFP_Data(ch).ntrode_id; fprintf('%i ',ch);
    tmpData = LFP_Data(ch).fields.data';
    
    % downsample the data
    if dsBool; tmpData = resampSig(tmpData, fs, fs_new); end
    
    lfpMat(ntrode, :) = tmpData;
    
    % for plotting the lfp sorting
    if showBuildMat
        if mod(ch,5) == 0
            subplot(1,6,subInd);
            imagesc(lfpMat(:,1:100));
            title(['Iter: ', num2str(ch)]);
            if ch == 5; ylabel('Sorting trodes LFP'); end
            subInd = subInd + 1;
        end
    end
end




%%  Wrap it all up with a bow and ribbon
lfpStruct = struct;
lfpStruct.data = lfpMat;
lfpStruct.info.dataDir =  dataDir; % todo change to meta
lfpStruct.info.path = path;
lfpStruct.info.dsBool = dsBool;
lfpStruct.info.fs = fs_new;
lfpStruct.info.date = date;
lfpStruct.info.origFormat = 'Trodes';
lfpStruct.info.low_pass_filter = LFP_Data(ch).low_pass_filter;
lfpStruct.info.voltage_scaling = LFP_Data(ch).voltage_scaling;
% lfpStruct.info.clockrate = LFP_Data(ch).clockrate; % TODO: NEED TO CONFIRM IF THIS IS THE FS!!!                                                                                                                                                                     
lfpStruct.info.original_file = LFP_Data(ch).original_file;

lfpStruct.info.ratName = metaRat.ratName; 

fprintf('\n');

%% Save our hard work!
% strip the .LFP/ from the file name
if saveStruct 
    structPath = [fullfile('D:','Documents','Data','lfpStructs'), filesep]; cd(structPath);
    save(lfpStruct.info.dataDir(19:end-5), 'lfpStruct'); 
    disp(['LFP Saved as ', structPath, '!']); 
end

fprintf('\nFinished! \n'); toc;

end