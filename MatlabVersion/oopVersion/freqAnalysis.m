function [freq,mag] = freqAnalysis(t,y)
% Take an FFT and return data

% t - time data , evenly spaced
% y - real space data

mag = abs(fft(y));
freq = linspace(0,1,length(y))/(t(2)-t(1));


end

