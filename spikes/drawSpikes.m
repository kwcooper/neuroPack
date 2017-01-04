function [] = drawSpikes(firingRate, duration, numNeurons)
[spikes, times] = makeSpikes(firingRate, duration, numNeurons);
plotSpikes(spikes, times*1000);
xlabel('Time (ms)');
ylabel('Neuron');
end

