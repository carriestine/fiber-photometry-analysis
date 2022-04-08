function [Zsignal] = getZscore(signal)

Zsignal = (signal-mean(signal))./std(signal);
