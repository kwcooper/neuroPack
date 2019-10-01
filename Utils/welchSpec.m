

function [f,logPower] = welchSpec(signal, fs, plt)
% Wrapper to compute the Welch PSD in Hz using matlab's built in f(x)
% docs: https://www.mathworks.com/help/signal/ref/pwelch.html
% ud 191001 KWC

rng default

% segment length of 500 samples with 300 overlapped samples
[pxx,f] = pwelch(signal,500,300,500,fs);
logPower = 10*log10(pxx);

if plt
    figure; plot(f,logPower);
    xlabel('Frequency (Hz)');
    ylabel('PSD (dB/Hz)');
end
end


