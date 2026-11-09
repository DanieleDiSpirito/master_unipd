function [lbp, nFloored] = logBandPower(sfilt, Rk, fs, winSec, powerFloor)
% LOGBANDPOWER Logarithmic band power: square -> moving average (FIR) -> ln.
%
%   [lbp, nFloored] = logBandPower(sfilt, Rk, fs)
%   [lbp, nFloored] = logBandPower(sfilt, Rk, fs, winSec, powerFloor)
%
%   sfilt      : [samples x channels] band-pass filtered EEG (e.g. bandpassRuns)
%   Rk         : [samples x 1] run index
%   fs         : sampling frequency in Hz
%   winSec     : moving average window in seconds (default 1)
%   powerFloor : positive value used in place of non-positive smoothed power
%                before the log (default eps). Should never be needed.
%
%   lbp        : [samples x channels] ln(P), P in uV^2 (natural log, not dB)
%   nFloored   : number of samples where the floor was applied
%
%   Square and moving average are done by bandPower (see there for the
%   transient at the start of each run).

    if nargin < 4 || isempty(winSec)
        winSec = 1;
    end
    if nargin < 5 || isempty(powerFloor)
        powerFloor = eps;
    end

    P = bandPower(sfilt, Rk, fs, winSec);

    nonPositive = P <= 0;
    nFloored = nnz(nonPositive);
    P(nonPositive) = powerFloor;

    lbp = log(P);

    if nFloored > 0
        warning('logBandPower: %d non-positive power samples set to %g', nFloored, powerFloor);
    end
    assert(all(isfinite(lbp(:))), 'Log band power has non-finite values');
end
