function [trials, Ck, RunTrial, info] = extractTrials(X, lbl, opts)
% EXTRACTTRIALS Cut the continuous signal(s) into trials using BOTH Tk and Rk.
%
%   [trials, Ck, RunTrial, info] = extractTrials(X, lbl)
%   [...] = extractTrials(X, lbl, StartFrom="fixation", MaxLengthDiff=2)
%
%   X    : [samples x channels] matrix, or cell array of matrices with the same
%          number of samples (e.g. {raw, logMu, logBeta}): every representation
%          gets exactly the same intervals and the same trimming
%   lbl  : label struct from labelEvents
%
%   Options
%       StartFrom     "trial" (default, trial start as in Tk) or "fixation"
%                     (fixation onset, as required from Lab04)
%       MaxLengthDiff max allowed difference (samples) between the longest
%                     and the shortest trial (default 2). Longer trials are
%                     trimmed at the END to the shortest length.
%
%   trials   : [samples x channels x trials] (cell of them if X is a cell)
%   Ck       : [trials x 1] cue of each trial
%   RunTrial : [trials x 1] run of each trial
%   info     : struct with, for each trial, start sample (in the concatenated
%              signal), original length, trial number inside the run, cue
%              and feedback onsets in samples relative to the trial start,
%              fixation / feedback = [first last] sample INDICES inside the
%              extracted trial ([trials x 2], 1-based, clipped by the trimming)

    arguments
        X
        lbl struct
        opts.StartFrom (1,1) string {mustBeMember(opts.StartFrom, ["trial", "fixation"])} = "trial"
        opts.MaxLengthDiff (1,1) double = 2
    end

    isCell = iscell(X);
    if ~isCell
        X = {X};
    end
    for k = 1:numel(X)
        assert(size(X{k}, 1) == numel(lbl.Tk), 'Signal %d and labels have different lengths', k);
    end

    runs = unique(lbl.Rk(lbl.Rk > 0))';

    starts = [];
    lens = [];
    Ck = [];
    RunTrial = [];
    TrialInRun = [];
    cueOnset = [];
    feedbackOnset = [];
    fixation = zeros(0, 2);
    feedback = zeros(0, 2);

    for r = runs
        trialNums = unique(lbl.Tk(lbl.Rk == r & lbl.Tk > 0))';
        for t = trialNums
            idx = find(lbl.Tk == t & lbl.Rk == r);
            assert(all(diff(idx) == 1), 'Trial %d of run %d is not contiguous', t, r);

            if opts.StartFrom == "fixation"
                fixIdx = find(lbl.Fk(idx) ~= 0, 1);
                assert(~isempty(fixIdx), 'Trial %d of run %d has no fixation', t, r);
                idx = idx(fixIdx:end);
            end

            startIdx = idx(1);
            endIdx = idx(end);
            assert(lbl.Rk(startIdx) == lbl.Rk(endIdx)); % interval inside the same run

            cues = unique(lbl.Ak(idx));
            cues(cues == 0) = [];
            assert(numel(cues) == 1, 'Trial %d of run %d has %d cues', t, r, numel(cues));

            starts(end+1, 1) = startIdx; %#ok<AGROW>
            lens(end+1, 1) = endIdx - startIdx + 1; %#ok<AGROW>
            Ck(end+1, 1) = cues; %#ok<AGROW>
            RunTrial(end+1, 1) = r; %#ok<AGROW>
            TrialInRun(end+1, 1) = t; %#ok<AGROW>
            cueOnset(end+1, 1) = find(lbl.Ak(idx) ~= 0, 1) - 1; %#ok<AGROW>
            feedbackOnset(end+1, 1) = find(lbl.CFk(idx) ~= 0, 1) - 1; %#ok<AGROW>

            fixRel = find(lbl.Fk(idx) ~= 0);
            fbRel = find(lbl.CFk(idx) ~= 0);
            assert(~isempty(fixRel) && all(diff(fixRel) == 1), 'Trial %d of run %d: bad fixation period', t, r);
            assert(~isempty(fbRel) && all(diff(fbRel) == 1), 'Trial %d of run %d: bad feedback period', t, r);
            fixation(end+1, :) = [fixRel(1) fixRel(end)]; %#ok<AGROW>
            feedback(end+1, :) = [fbRel(1) fbRel(end)]; %#ok<AGROW>
        end
    end

    minLen = min(lens);
    maxLen = max(lens);
    if maxLen - minLen > opts.MaxLengthDiff
        error(['Trial lengths go from %d to %d samples (diff > %d): ' ...
            'check the extraction before trimming'], minLen, maxLen, opts.MaxLengthDiff);
    end

    % periods are cut by the trimming at the END too
    fixation = min(fixation, minLen);
    feedback = min(feedback, minLen);
    assert(all(fixation(:, 1) <= fixation(:, 2)) && all(feedback(:, 1) <= feedback(:, 2)), ...
        'A period is empty after trimming');

    nTrials = numel(starts);
    trials = cell(size(X));
    for k = 1:numel(X)
        trials{k} = zeros(minLen, size(X{k}, 2), nTrials);
        for j = 1:nTrials
            trials{k}(:, :, j) = X{k}(starts(j):starts(j) + minLen - 1, :); % trim at the END
        end
    end
    if ~isCell
        trials = trials{1};
    end

    info = struct('start', starts, 'length', lens, 'trialInRun', TrialInRun, ...
        'cueOnset', cueOnset, 'feedbackOnset', feedbackOnset, ...
        'fixation', fixation, 'feedback', feedback, ...
        'minLength', minLen, 'maxLength', maxLen);
end
