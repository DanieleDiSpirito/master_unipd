function values = periodMean(trials, period)
% PERIODMEAN Mean of each trial over its own period (e.g. feedback or fixation).
%
%   values = periodMean(trials, period)
%
%   trials : [samples x channels x trials]
%   period : [trials x 2] [first last] sample indices inside each trial
%            (e.g. info.feedback / info.fixation from extractTrials)
%
%   values : [trials x channels] one value per channel and trial

    nTrials = size(trials, 3);
    assert(size(period, 1) == nTrials, 'period and trials describe a different number of trials');

    values = zeros(nTrials, size(trials, 2));
    for j = 1:nTrials
        values(j, :) = mean(trials(period(j, 1):period(j, 2), :, j), 1);
    end
end
