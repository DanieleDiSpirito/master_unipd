function [ERD_ERS, Reference, valid] = erdErs(TrialPower, refPeriod)
% ERDERS Percentage ERD/ERS of each trial w.r.t. its own reference period.
%
%   [ERD_ERS, Reference, valid] = erdErs(TrialPower, refPeriod)
%
%   TrialPower : [samples x channels x trials] smoothed power (NOT log)
%   refPeriod  : [trials x 2] [first last] reference samples inside each
%                trial (e.g. info.fixation from extractTrials)
%
%   ERD_ERS    : [samples x channels x trials] 100 * (P - R) ./ R
%                (negative = ERD, positive = ERS). NaN for invalid trials.
%   Reference  : [1 x channels x trials] mean power in the reference period
%   valid      : [trials x 1] true if every reference value of the trial is
%                finite and strictly positive

    nTrials = size(TrialPower, 3);
    assert(size(refPeriod, 1) == nTrials, 'refPeriod and TrialPower describe a different number of trials');

    Reference = zeros(1, size(TrialPower, 2), nTrials);
    for j = 1:nTrials
        Reference(1, :, j) = mean(TrialPower(refPeriod(j, 1):refPeriod(j, 2), :, j), 1);
    end

    valid = squeeze(all(isfinite(Reference) & Reference > 0, 2));
    valid = valid(:);

    ERD_ERS = 100 * (TrialPower - Reference) ./ Reference; % implicit expansion along samples
    ERD_ERS(:, :, ~valid) = NaN;
end
