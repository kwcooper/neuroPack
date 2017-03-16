function [spikeMatDist, tVec] = distroSpikeGen(tSim, nTrials, fl, fh)
%Spike generator
%input time steps(in s), number of trials, and lower and upper frequencys
%output an array of spikes and time Vectors

dt = 1/1000; % convert to ms steps
nBins = floor(tSim/dt);
spikeMatDist = zeros(nTrials,nBins);

%create a numrical distribution 
%dist = [10:nBins/2:30 fliplr(10:(nBins/2):30)];
dist = [linspace(fl, fh, nBins/2) fliplr(linspace(fl, fh, nBins/2))];

ii = 1;
for fr = dist
    spikeMatDist(:,ii) = rand(nTrials, 1) < fr*dt;
    ii = ii + 1;
end
tVec = 0:dt:tSim-dt;

spikeMatDist = logical(spikeMatDist);
end




