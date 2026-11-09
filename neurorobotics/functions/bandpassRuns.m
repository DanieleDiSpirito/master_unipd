function [sfilt, b, a] = bandpassRuns(s, Rk, band, order, fs)
% BANDPASSRUNS Zero-phase Butterworth band-pass filter applied SEPARATELY to each run.
%
%   [sfilt, b, a] = bandpassRuns(s, Rk, band, order, fs)
%
%   s     : [samples x channels] continuous EEG
%   Rk    : [samples x 1] run index (from labelEvents)
%   band  : [low high] cutoff frequencies in Hz (e.g. [8 13])
%   order : Butterworth order passed to butter (band-pass -> 2*order poles)
%   fs    : sampling frequency in Hz
%
%   sfilt : [samples x channels] filtered EEG, same positions as s
%   b, a  : filter coefficients (e.g. for filterAnalyzer(b, a, 'SampleRates', fs))
%
%   The join between two runs is not real continuous EEG: filtering across
%   it would introduce artificial transients, so each run is filtered alone.

    assert(numel(band) == 2 && band(1) > 0 && band(1) < band(2) && band(2) < fs / 2, ...
        'Cutoffs [%g %g] Hz must satisfy 0 < low < high < fs/2 = %g Hz', band(1), band(2), fs / 2);
    assert(size(s, 1) == numel(Rk), 'Signal and Rk have different lengths');

    [b, a] = butter(order, 2 * band / fs);

    sfilt = zeros(size(s));
    for r = unique(Rk(Rk > 0))'
        idx = Rk == r;
        sfilt(idx, :) = filtfilt(b, a, s(idx, :)); % filtfilt works column-wise
    end

    assert(isequal(size(sfilt), size(s)));
    assert(all(isfinite(sfilt(:))), 'Filtered signal has non-finite values');
end
