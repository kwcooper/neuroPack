
%figPath = fullfile(dropboxPath,SLpath,'..','figures');

%this is keiland's version

% update these values or copy and paste them
Rat = 'Tio';
Session = '170703_1251_CircleTrack';
Recording = '2017-07-03_13-08-30';
tioChOrd = [43 44 46 45 40 39 37 38 59 60 58 57 52 51 53 54]; 
chOrdTxt = 'Probe order';
chOrd = tioChOrd;

% Rat = 'Rio';
% Session = '2017-08-10_CircleTrack';
% Recording = '2017-08-10_19-14-01';
% rioChOrd = [43 44 46 45 40 39 37 38 59 60 58 57 52 51 53 54]; 
% chOrdTxt = 'Probe order';
%chOrd = rioChOrd;

% update this to match the name of the variable with your channel order

dsFreq = 600; %data sample frequency? what is this?
ref = 8; %channel reference
nChan = length(chOrd);
fsVid = 120;

workingDir = fullfile(ratLibPath,Rat,Session,Recording); cd(workingDir);


%% 
if ~exist('root', 'var')
    
    %[data, timestamps, info] = load_open_ephys_data_faster(filename, varargin)
    
    tmpRoot = CMBHOME.Session('name','experiment1','epoch',[-inf inf],'b_ts',[0:0.01:21.0028*60],'fs_video',120);
    tmpRoot.path_lfp = repmat({'experiment1_100.raw.kwd'},1,1);
    tmpRoot = tmpRoot.LoadLFP(1,'downsample',dsFreq,'chOrd',chOrd(1));
    tmpRoot.active_lfp = nChan;
    
    
    fprintf('Creating root object... \n')
    root = CMBHOME.Session('name','experiment1','epoch',[-inf inf],'b_ts',[0:1/fsVid:max(tmpRoot.b_lfp(1).ts)],'fs_video',fsVid);
    root.path_lfp = repmat({'experiment1_100.raw.kwd'},1,nChan);
    %save experiment1.mat root
    
    fprintf('Loading lfp... \n')
    root = root.LoadLFP(1:nChan,'downsample',dsFreq,'chOrd',chOrd);
    root.active_lfp = nChan;
end

% %if isempty(root.b_lfp(1)) %need to fix this somehow
%     fprintf('Loading lfp... \n')
%     root = root.LoadLFP(1:nChan,'downsample',dsFreq,'chOrd',chOrd);
%     root.active_lfp = nChan;
% %end

root = rmDataBlips(root);

%iterate through the channels
chans = [1:nChan];
dataDS = nan(length(chans),length(CMBHOME.Utils.ContinuizeEpochs(root.lfp.signal)));
for i = 1:length(chans)
  root.active_lfp = chans(i);
  dataDS(i,:) = CMBHOME.Utils.ContinuizeEpochs(root.lfp.signal);
end

%% Extract Theta

Fs = root.lfp.fs;
tInfo = {};
tInfo.Wn_theta = [6/(Fs/2) 10/(Fs/2)];
[tInfo.btheta,tInfo.atheta] = butter(3,tInfo.Wn_theta);


% extract theta with phase and power
fprintf('theta extraction \n')

tInfo.theta_filt = nan(size(dataDS));
tInfo.theta_phase =  nan(size(dataDS));
tInfo.theta_amp =  nan(size(dataDS));
% why is this iterating through each data point, and not taking all of the data? k
for iD =  1:size(dataDS,1)
  tInfo.theta_filt(iD,:) = filtfilt(tInfo.btheta,tInfo.atheta,dataDS(iD,:)); %filter the data
  tInfo.theta_phase(iD,:) = atan2(imag(hilbert(tInfo.theta_filt(iD,:))), tInfo.theta_filt(iD,:));
  tInfo.theta_amp(iD,:) = abs(hilbert(tInfo.theta_filt(iD,:)));
end

%% 
% Focus on epochs of data with high theta amplitude
% high amplitude will be based on greater than 2 std above the mean
% so, compute session-wide mean and standard deviation of theta power
% using the last channel (update this if another channel makes more sense)

% I might need to change this k 
% Also, why not all theta? quick n dirty or the whole project?
ch = size(tInfo.theta_amp,1);
tInfo.thetaMeanAmp = mean(tInfo.theta_amp(ch,:));
stdAmp = std(tInfo.theta_amp(ch,:));
 
% now find the high theta 
%highTheta = find(theta_amp(ch,:)>(meanAmp+2*stdAmp));
%highTheta = find(theta_amp(end,:)>(meanAmp));
tInfo.highTheta = find(tInfo.theta_amp(end,:)>(tInfo.thetaMeanAmp-stdAmp));
tInfo.highThetaEp = mat2cell(tInfo.highTheta, 1, diff([0 find([(diff(tInfo.highTheta) > 1) 1])]));
lengthEp = cellfun(@length,tInfo.highThetaEp);
inds_long = lengthEp>300;
tInfo.highThetaEp_long = tInfo.highThetaEp(inds_long);
 
%find the theta shift for high theta channels
clear thetaShift
for iT = 1:length(tInfo.highThetaEp_long)
  for iC1 = 1:length(chans)
    for iC2 = 1:length(chans)
      tInfo.thetaShift{iT}(iC1,iC2,:) = circDiff([tInfo.theta_phase(iC1,tInfo.highThetaEp_long{iT})', ...
                                                  tInfo.theta_phase(iC2,tInfo.highThetaEp_long{iT})'],2,'rad'); % changed orientation of data
    end
  end
end

tInfo.thetaShiftMat = cat(3,tInfo.thetaShift{:}); % used cat instead of cell2mat
fprintf('Basing calculations off of %2.2f s of data\n',size(tInfo.thetaShiftMat,3)/Fs);
[tInfo.thetaShiftAngle,tInfo.thetaShiftRbar] = circmean(tInfo.thetaShiftMat,3);


%%
%new theta Extraction
disp("Extracting theta cycles...");
metho = 'hilbert';
disp("useing " + metho)
root.epoch=[-inf,inf];
band = [6,10];
[thetaPhs,~,~] = extractThetaPhase(root.b_lfp(16).signal,Fs,metho,band);
[cycles,~] = parseThetaCycles(thetaPhs,Fs,band);

inds = find(cycles);
%%
%Figures
%keyboard
quiverPlot(tInfo,chOrdTxt)
corrPlot(tInfo,chOrdTxt)
thetaGreaterMeanPower(tInfo, chans)
%plotAvgWaveImg(root, cycles, ref)
plotAvgWaves(root,Fs)
plotRawWaves(root,Fs)



function quiverPlot(tInfo,chOrdTxt)
% plots the polar coherence between electrodes
[u,v] = pol2cart(tInfo.thetaShiftAngle,tInfo.thetaShiftRbar);
figure; quiver(u(2:2:end,2:2:end),v(2:2:end,2:2:end)); axis ij;  
title([chOrdTxt ', Blank (1403)']);
end

function corrPlot(tInfo, chOrdTxt)
%plots the coherence between each cell
R = corr(tInfo.theta_filt');
figure; imagesc(R);
FirstDegMCorr = mean(diag(R,1));
SecondDegMCorr = mean(diag(R,2)); 
title([chOrdTxt, ', Blank (1403), 1st vs 2nd neighbor corr ratio = ',num2str(FirstDegMCorr/SecondDegMCorr)]);

[u,v] = pol2cart(tInfo.thetaShiftAngle,tInfo.thetaShiftRbar);
circR = sqrt(u.^2 + v.^2); figure; imagesc(circR,[0 1]);
title("corr in radians");
end

function thetaGreaterMeanPower(tInfo, chans)
%highTheta = find(theta_amp(end,:)>(meanAmp+2*stdAmp));
tInfo.highTheta = find(tInfo.theta_amp(end,:)>(tInfo.thetaMeanAmp));
tInfo.highThetaEp = mat2cell(tInfo.highTheta, 1, diff([0 find([(diff(tInfo.highTheta) > 1) 1])]));
lengthEp = cellfun(@length,tInfo.highThetaEp);
inds_long = lengthEp>300;
tInfo.highThetaEp_long = tInfo.highThetaEp(inds_long);
 
clear thetaShift
for iT = 1:length(tInfo.highThetaEp_long),
  tInfo.thetaShift{iT} = circDiff(tInfo.theta_phase(:,tInfo.highThetaEp_long{iT}),1,'rad');
end
tInfo.thetaShiftMat = cell2mat(tInfo.thetaShift);
 
bins = [-pi:0.05:pi];
figure; 
for i = 1:size(tInfo.thetaShiftMat,1), 
  subplot(size(tInfo.thetaShiftMat,1),1,i), hist(tInfo.thetaShiftMat(i,:),bins); 
  if i==1, title('8/17/2016 14:03 Recordings. Theta > mean power.'); end
  ylabel([num2str(chans(i+1)) ' - ' num2str(chans(i))]);
end
end

function plotAvgWaveImg(root, cycles, ref)
CycleTs=root.b_lfp(ref).ts(cycles);
Epochs = [CycleTs-0.125 CycleTs+0.125];
root.epoch=Epochs;

%for each channel, fetch the eeg snips, then average them
for I=1:16
    root.active_lfp=I; 
    EegSnips=root.lfp.signal;
    MeanThetaWave(I,:)=nanmean(catPlus(3,EegSnips),3);
end

figure; imagesc(MeanThetaWave)
title("Mean Theta Wave; ref:" + string(ref))

end

function plotRawWaves(root, Fs)
%make the raw waves from lfp signal
for i=1:16
    rawWaves(i, :) = root.b_lfp(i).signal;
end
%plot raw eegWaves
kPlotLFP(rawWaves(8:2:16,:),1:1200,Fs)

%plotLFP(rawWaves(8:2:16,:), Fs);

end

function plotAvgWaves(root, Fs)
for I=1:16
    root.active_lfp=I;
    EegSnips=root.lfp.signal;
    MeanThetaWave(I,:)=nanmean(catPlus(3,EegSnips),3);
end
data = MeanThetaWave(8:2:end,:);
%plotLFP(data, Fs, 2.5)
kPlotLFP(data,1:1200,Fs, 2.5)
end
