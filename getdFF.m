function [dFFsignal,dFFbaselinesignal] = getdFF(signal, reference, basepersignal)

% Calculates DF/F ((F-F0)/F0) for a signal using the median of the entire trace as
% F0, and/or signal using median of a defined baseline period as F0.


if nargin < 3
    basepersignal = signal;
    if nargin < 2
        reference = signal;
    end
end


dFFsignal = (signal - (median(reference)))./abs(median(reference));
%dFFsignal = dFFsignal.*100;

dFFbaselinesignal = (signal-(median(basepersignal)))./abs(median(basepersignal));
dFFbaselinesignal = dFFbaselinesignal.*100; %turn into percentage