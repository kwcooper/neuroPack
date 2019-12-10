
function plx2nex()
%% A script to load .plx filex and convert  them into nex files
% we use this mainly to deal with the alignment issues in opensorter. 
% requires a handful of libs and files to properly run

% Look towards convertPLXtoMDA.m and ConvertAnalyzedMDAtoNEX.m for more info
% https://github.com/FortinLab/MSssInteg

% TODO:
%   Need to add logic for file-wise vs tetrode-wise

% init KWC 191207; last updated KWC 191209

%%
load('RosettaLC.mat'); % used for naming channels
% Params for the waveforms -> derived from the sample JSON file
params.samplerate = 40000;
params.clip_size = 32; 
tetName = 'Stella_novel_1';


%% Identify PLX file
[fileName, filePath] = uigetfile('.plx','Identify .PLX File');
if fileName == 0
    disp('No file selected');
    return
end
plxFile = [filePath fileName];
origDir = cd; cd(filePath);

mkdir(tetName);  % create a new directory to store the files
fDir = tetName;

%% Identify channels and waveforms
% Pull out relevant information from the recording file
%   sampRate = sample rate the recording was done at
%   nTrodes = number of channels per recording, tetrode = 4
%   preThresh = number of samples recorded prior to threshold crossing for
%       trace capture
%   ssnDur = duration of the recording in seconds
[~,~,sampRate,~,nTrodes,~,preThresh,~,~,~,~,ssnDur,~] = plx_information(plxFile);


[tsCountFl, ~, ~, ~] = plx_info(plxFile, 1);
tetChans = 1:nTrodes:(size(tsCountFl,2));

nexFile = nexCreateFileData(params.samplerate); 
trode = 1;
wireCnt = 1;
for chan = 1:size(tsCountFl,2)-1
    tic
    % Pull out waveform data from the recording. NOTE: The input of 0 means
    % pull the unsorted units so make sure you're only using an UNCUT file
    % in this analysis. OR that you've put all the crappy waveforms into
    % units and have the CLEAN data in the UNSORTED (0) unit.
    %   numWFs = number of waveforms on that channel
    %   npw = number of points recorded per wave
    %   ts = timestamps associated with each waveform
    %   wave = waveform values in mV
    [numWFs, npw, ts, wave] = plx_waves_v(plxFile, chan, 0);
    waves(chan) = {wave};
    
    % check for empty waves
    if isequal(wave, -1)
        wave = [];
        ts = [];
    end
    
    if 1
        % nexFile, WFreq, timestamps, waveforms, 
        % name, preThresholdTimeInSeconds, 
        % numberOfPointsInWaveform, wireNumber, unitNumber 
    nexFile = nexAddWaveform(nexFile, params.samplerate, ts, wave',...
                sprintf('%s_%s_wf_%d', tetName, RosettaLC{trode}, wireCnt), preThresh/params.samplerate,...
                params.clip_size, chan, 0);
    end
    
    % Save each tetrode to it's own file
    fName = [tetName,'_',RosettaLC{trode},'.nex'];
    fPath = fullfile(fDir, fName);
    res = writeNexFile(nexFile, fPath); 
    if res
        fprintf('%s.nex Saved\n',tetName);
    else
        disp('Uh Oh...')
    end
    
    % Add the logic for figuring out which of the tetrodes and wires we're on
    if isequal(wireCnt, 4)
        wireCnt = 1;
        trode = trode + 1;
        nexFile = nexCreateFileData(params.samplerate); 
        disp(['on trode' num2str(trode)])
    else
        wireCnt = wireCnt + 1;
    end
    
    toc
end
end





