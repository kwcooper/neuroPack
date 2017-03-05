function [spikes, times] = makeSpikes(firingRate, duration, numNeurons)
%Dayan & Abbot 30
%Thanks to RAJAT SAXENA

dt = 1/1000; % s
chunks = floor(duration/dt);
spikes = rand(numNeurons, chunks) < firingRate*dt;
times = 0:dt:duration-dt;