   
function [freq,mag] = freqMag(x, Fs) 
% Simple func to grab the frequency spectrum for a signal (fft)
%
% INPUTS: 
% x = signal vector
% Fs = sampleing rate
% 
% EXAMPLE CALL:
% Grab the spectrum for a channel
% [freq,mag] = freqMag(lfpStruct.data(1,:), fs); 
% 
% %%Plot the spectrum:
% figure;
% plot(freq,mag);
% xlabel('Frequency (Hz)');
% title('Magnitude Response');
% xlim([-1 90]); % set the limits for the frequencies you want to use


% Time specifications:
dt = 1/Fs;                     % seconds per sample
StopTime = length(x)/Fs;                  % seconds
t = (0:dt:StopTime-dt)';
N = size(t,1);

% Fourier Transform:
X = fftshift(fft(x));

% Frequency specifications:
dF = Fs/N;                      % hertz
freq = -Fs/2:dF:Fs/2-dF;           % hertz

mag = abs(X)/N;



