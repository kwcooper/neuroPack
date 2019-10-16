function xNew = resampSig(x, Fs1, Fs2)
% resample a signal using Matlab's resample function
% calculates the integer values needed from the signal's
% sample rate (Fs1) to the desired sample rate (Fs2). 

[P,Q] = rat(Fs2/Fs1);
xNew = resample(x,P,Q);

if 0 % validate that the signals and time is conserved (disabled)
t1 = linspace(0,size(x,2)/Fs1,size(x,2));
figure; plot(t1,x)
t2 = linspace(0,size(xNew,2)/Fs2,size(xNew,2));
figure; plot(t2,xNew)
end

end