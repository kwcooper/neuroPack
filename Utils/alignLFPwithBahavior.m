function alignLFPwithBahavior(behavData, lfpStruct, buff_secs)
%
%
%  The timestamps are computed by:
%    curTime = clock;
%    %            Hours               Minutes          Seconds
%    timeStamp = (curTime(4)*3600) + (curTime(5)*60) + curTime(6);
%
%   A la Gabes' function
%
%
%

viz = 0; % to see the new alignment

% Grab the behav data and the lfp
bd = behavData.ssnData;
lfp = lfpStruct.data(10,:); % TODO change

buff_secs = 122.13; % add buffer to task data in secs
fs = lfpStruct.info.fs; 

% compute some LFP stats that we may want later
LFP_length = size(lfp,2);
LFP_secs = LFP_length / fs;
LFP_min = LFP_secs / 60;

timeax = linspace(0,LFP_secs,LFP_length); % compute the lfp time axis


% collect the poke in and out times
pits = [bd.PokeInTime]; pots = [bd.PokeOutTime];

% Right now we just want it algigned via seconds. So we round. 
% in the future if we need more accuracy then we will have to
% be more specific about how we compute this. 
% also normalize into seconds from zero, not the arbitary 
% task start time
rpits_redux = round(pits - pits(1));
rpots_redux = round(pots - pots(1));


task_min = (rpits(end) - rpits(1) ) / 60; % This isn't the most accurate, use poke out
task_fs = 1; % TODO

% need to generate the poking data at the same fs as the lfp
tstPke = linspace(0, rpits_redux(end), rpits_redux(end)*fs);
pkes = ismember(round(tstPke), rpits_redux);
buf = zeros(1,buff_secs*fs);
newData = [buf, pkes];  % add the offset to the poking data
pkeTme = linspace(0, size(pkes,2)/fs, size(pkes,2));
bufTme = linspace((-size(buf,2)/fs), 0, size(buf,2)+1);
behavTme = [bufTme(:,1:end-1), pkeTme];

% Check if what you get is what you want
if viz
    figure; 
    subplot(3,1,1);
    plot(pkeTme)
    subplot(3,1,2);
    plot(bufTme)
    subplot(3,1,3);
    plot(behavTme)
end

% now let's look at the 

end