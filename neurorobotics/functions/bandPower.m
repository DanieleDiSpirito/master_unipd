function [P, sq] = bandPower(sfilt, Rk, fs, winSec)
% BANDPOWER Smoothed band power: square -> moving average (FIR, filter), NO log.
%
%   [P, sq] = bandPower(sfilt, Rk, fs)
%   [P, sq] = bandPower(sfilt, Rk, fs, winSec)
%
%   sfilt  : [samples x channels] band-pass filtered EEG (e.g. bandpassRuns)
%   Rk     : [samples x 1] run index
%   fs     : sampling frequency in Hz
%   winSec : moving average window in seconds (default 1)
%
%   P      : [samples x channels] smoothed power in uV^2 (input of ERD/ERS)
%   sq     : [samples x channels] squared EEG (instantaneous power)
%
%   The moving average is a trailing window (current + previous samples)
%   computed with filter, restarted at every run: the first winSec seconds
%   of each run contain the initialization transient and the output is
%   delayed by about winSec/2.

    if nargin < 4 || isempty(winSec)
        winSec = 1;
    end
    assert(size(sfilt, 1) == numel(Rk), 'Signal and Rk have different lengths');

    windowSize = round(winSec * fs);
    b = ones(1, windowSize) / windowSize; % equal coefficients, sum = 1
    a = 1;

    sq = sfilt .^ 2; % instantaneous power
    P = zeros(size(sfilt));
    for r = unique(Rk(Rk > 0))'
        idx = Rk == r;
        P(idx, :) = filter(b, a, sq(idx, :), [], 1); % new call -> filter state reset for this run
    end

    assert(isequal(size(P), size(sfilt)));
    assert(all(isfinite(P(:))), 'Band power has non-finite values');
    assert(all(P(:) >= 0), 'Band power has negative values');
end
