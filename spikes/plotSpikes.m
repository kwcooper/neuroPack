function [] = plotSpikes(spikes, times)
hold all;
for iterations = 1:size(spikes,1)
    place = times(spikes(iterations, :));
    for count = 1:length(place)
        plot([place(count) place(count)], [iterations-0.4 iterations+0.4], 'k');
    end
end
ylim([0 size(spikes, 1)+1]);