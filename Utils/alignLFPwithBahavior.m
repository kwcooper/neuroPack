function [alignedAx] = alignLFPwithBahavior(behavData, lfpStruct, buff_secs)
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

%buff_secs = 122.13; % add buffer to task data in secs
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
buf = zeros(1,buff_secs*fs);

% Shift the poke in/out data according to the buffer
% need to generate the poking data at the same fs as the lfp
pits_empty = linspace(0, rpits_redux(end), rpits_redux(end)*fs);
pits_log = ismember(round(pits_empty), rpits_redux);
pits_new = [buf, pits_log];  % add the offset to the poking data

pots_empty = linspace(0, rpots_redux(end), rpots_redux(end)*fs);
pots_log = ismember(round(pots_empty), rpots_redux);
pots_new = [buf, pots_log];  % add the offset to the poking data

% compute the poke in/out time axis
bufTme = linspace((-size(buf,2)/fs), 0, size(buf,2)+1);
pits_tme = linspace(0, size(pits_log,2)/fs, size(pits_log,2));
pitsTmeax = [bufTme(:,1:end-1), pits_tme];
pots_tme = linspace(0, size(pots_log,2)/fs, size(pots_log,2));
potsTmeax = [bufTme(:,1:end-1), pots_tme];

% now let's fix the lfp time axis
lfpTmeax = [bufTme(:,1:end-1), timeax];
lfpTmeax = lfpTmeax(:,1:size(lfp,2));

% Check if what you get is what you want
if viz
    figure; 
    a = subplot(3,1,1);
    plot(lfpTmeax, lfp)

    b = subplot(3,1,2);
    plot(pitsTmeax, pits_new)
    
    c = subplot(3,1,3);
    plot(potsTmeax, pots_new)
    linkaxes([a,b,c],'x');
end

% pack it all up
alignedAx.lfpTmeax = lfpTmeax;
alignedAx.lfp = lfp;
alignedAx.pitsTmeax = pitsTmeax;
alignedAx.pits_new = pits_new;
alignedAx.potsTmeax = potsTmeax;
alignedAx.pots_new = pots_new;
alignedAx.stats.task_min = task_min;
alignedAx.stats.task_fs = task_fs;

end