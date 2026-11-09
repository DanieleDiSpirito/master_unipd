function [avg, counts, cues, sem] = cueAverage(trials, Ck, cues)
% CUEAVERAGE Average (and standard error) of the trials of each cue (across all runs).
%
%   [avg, counts, cues, sem] = cueAverage(trials, Ck)
%   [avg, counts, cues, sem] = cueAverage(trials, Ck, [771 773])
%
%   trials : [samples x channels x trials]
%   Ck     : [trials x 1] cue of each trial
%   cues   : cues to average (default unique(Ck))
%
%   avg    : [samples x channels x numel(cues)]
%   counts : [numel(cues) x 1] trials contributing to each average
%   sem    : [samples x channels x numel(cues)] standard error of the mean
%            = sample std across trials / sqrt(counts)

    if nargin < 3 || isempty(cues)
        cues = unique(Ck);
    end
    assert(size(trials, 3) == numel(Ck), 'trials and Ck describe a different number of trials');

    avg = zeros(size(trials, 1), size(trials, 2), numel(cues));
    sem = zeros(size(avg));
    counts = zeros(numel(cues), 1);
    for i = 1:numel(cues)
        sel = Ck == cues(i);
        counts(i) = nnz(sel);
        avg(:, :, i) = mean(trials(:, :, sel), 3);
        sem(:, :, i) = std(trials(:, :, sel), 0, 3) / sqrt(counts(i)); % 0 -> sample std (N-1)
    end
end
