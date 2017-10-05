
%Plots the raw data under the extracted cycle
% figure; plot(root.b_lfp(1).signal(1:600));
% hold on; plot(100*thetaPhs(1:600));
% X=find(cycles(1:600)); plot(X,thetaPhs(X),'rx');
% title("waveform")

%%
metho = 'hilbert';
disp("useing" + string(metho));
Fs = 600;
ref = 8;


root.epoch=[-inf,inf];
band = [6,10];
[thetaPhs,~,~] = extractThetaPhase(root.b_lfp(ref).signal,Fs,metho,band);
[cycles,~] = parseThetaCycles(thetaPhs,Fs,band);

%%
%Grab the cycles in lfp
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

% figure;
% for I=1:16
%     hold on; plot(MeanThetaWave(I,:));
%     %find(cycles)?
%     legendInfo{I} = [num2str(I)];
% end
% title("Sacked waves" + string(ref))
% legend(legendInfo);

%make the raw waves from lfp signal
for i=1:16
    rawWaves(i, :) = root.b_lfp(i).signal;
end

%plot raw eegWaves
plotLFP(rawWaves(8:2:16,:), Fs);



