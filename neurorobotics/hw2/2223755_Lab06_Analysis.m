clear; close all; clc;

%% Configuration

MAT_FILES = { ...
    'ah7.20170613.161402.offline.mi.mi_bhbf.mat', ...
    'ah7.20170613.162331.offline.mi.mi_bhbf.mat', ...
    'ah7.20170613.162934.offline.mi.mi_bhbf.mat'};
OUT_PDF = '2223755_Lab06.pdf';

CHANNELS = {'C3', 'Cz', 'C4'};
cues = [773 771];
cueNames = {'both hands', 'both feet'};
MAX_LENGTH_DIFF = 2;

%% 4.1

N_RUNS = numel(MAT_FILES);
runs = cell(N_RUNS, 1);
for i = 1:N_RUNS
    runs{i} = load(MAT_FILES{i});
end

params = {'f', 'chLabels', 'samplerate', 'wlength', 'wshift', 'pshift', 'mlength', 'winconv'};
for i = 2:N_RUNS
    for k = 1:numel(params)
        assert(isequal(runs{i}.(params{k}), runs{1}.(params{k})), 'Run %d: different %s', i, params{k});
    end
    assert(size(runs{i}.PSD, 2) == size(runs{1}.PSD, 2) && size(runs{i}.PSD, 3) == size(runs{1}.PSD, 3));
end

f = runs{1}.f;
chLabels = runs{1}.chLabels;
wshift = runs{1}.wshift;
fprintf('Spectrogram: wlength %g s, wshift %g s, pshift %g s, mlength %g s, conversion ''%s''\n', ...
    runs{1}.wlength, wshift, runs{1}.pshift, runs{1}.mlength, runs{1}.winconv);
fprintf('Frequencies: %g..%g Hz (%d bins)\n', f(1), f(end), numel(f));

%% 4.2

sizesWin = cellfun(@(x) size(x.PSD, 1), runs);
offsetsWin = [0; cumsum(sizesWin(1:end-1))];

PSD = [];
POS = [];
DUR = [];
TYP = [];
for i = 1:N_RUNS
    PSD = cat(1, PSD, runs{i}.PSD);
    POS = [POS; runs{i}.EVENTwin.POS + offsetsWin(i)];
    DUR = [DUR; runs{i}.EVENTwin.DUR];
    TYP = [TYP; runs{i}.EVENTwin.TYP];
end
EVENTwin = struct('POS', POS, 'DUR', DUR, 'TYP', TYP);

lbl = labelEvents(EVENTwin, sizesWin);

%% 4.3 - 4.4

[nWinTot, nFreq, nCh] = size(PSD);
[trials, Ck, RunTrial, info] = extractTrials(reshape(PSD, nWinTot, nFreq * nCh), lbl, MAX_LENGTH_DIFF);
nWin = size(trials, 1);
nTrials = numel(Ck);
Activity = reshape(trials, nWin, nFreq, nCh, nTrials);

assert(size(unique([RunTrial info.trialInRun], 'rows'), 1) == nTrials);
fprintf('Trials: %d, length %d..%d windows -> %d (%.2f s)\n', nTrials, info.minLength, info.maxLength, nWin, nWin * wshift);
for i = 1:N_RUNS
    fprintf(' - run %d: %d hands, %d feet\n', i, nnz(RunTrial == i & Ck == 773), nnz(RunTrial == i & Ck == 771));
end

refPeriod = info.fixation;
actPeriod = info.feedback;
fprintf('Reference = fixation, windows %d..%d to %d..%d\n', min(refPeriod(:, 1)), max(refPeriod(:, 1)), min(refPeriod(:, 2)), max(refPeriod(:, 2)));
fprintf('Activity = continuous feedback, windows %d..%d to %d..%d\n', min(actPeriod(:, 1)), max(actPeriod(:, 1)), min(actPeriod(:, 2)), max(actPeriod(:, 2)));

%% 5

R = zeros(1, nFreq, nCh, nTrials);
for j = 1:nTrials
    R(1, :, :, j) = mean(Activity(refPeriod(j, 1):refPeriod(j, 2), :, :, j), 1);
end

valid = squeeze(all(all(isfinite(R) & R > 0, 2), 3));
fprintf('Excluded trials: %d\n', nnz(~valid));
for j = find(~valid)'
    fprintf(' - run %d, trial %d\n', RunTrial(j), info.trialInRun(j));
end

ERD_ERS = 100 * (Activity - R) ./ R;
ERD_ERS = ERD_ERS(:, :, :, valid);
Ck = Ck(valid);
cueOn = info.cueOnset(valid);
fbOn = info.feedbackOnset(valid);

%% 6

chIdx = zeros(1, numel(CHANNELS));
for c = 1:numel(CHANNELS)
    chIdx(c) = find(strcmpi(chLabels, CHANNELS{c}));
end

maps = zeros(nWin, nFreq, numel(cues), numel(CHANNELS));
counts = zeros(numel(cues), 1);
for k = 1:numel(cues)
    counts(k) = nnz(Ck == cues(k));
    for c = 1:numel(CHANNELS)
        maps(:, :, k, c) = mean(ERD_ERS(:, :, chIdx(c), Ck == cues(k)), 4);
    end
end

t = (0:nWin-1) * wshift;
fprintf('Cue onset %.3f..%.3f s, feedback onset %.3f..%.3f s from fixation\n', ...
    min(cueOn) * wshift, max(cueOn) * wshift, min(fbOn) * wshift, max(fbOn) * wshift);

L = prctile(abs(maps(:)), 99);
fprintf('Color limits: +-%.0f%% (99th percentile, max %.0f%%)\n', L, max(abs(maps(:))));

x = linspace(0, 1, 128)';
cmap = [x x ones(128, 1); ones(128, 1) flipud(x) flipud(x)];

fig = figure('Color', 'w', 'Position', [50 50 1500 800]);
for k = 1:numel(cues)
    for c = 1:numel(CHANNELS)
        ax = subplot(numel(cues), numel(CHANNELS), (k - 1) * numel(CHANNELS) + c);
        imagesc(ax, t, f, maps(:, :, k, c)');
        axis(ax, 'xy');
        colormap(ax, cmap);
        clim(ax, [-L L]);
        cb = colorbar(ax);
        cb.Label.String = 'ERD/ERS (%)';
        hold(ax, 'on');
        xline(ax, median(cueOn) * wshift, 'k--', 'Cue');
        xline(ax, median(fbOn) * wshift, 'k--', 'Feedback');
        hold(ax, 'off');
        title(ax, sprintf('%s - %d %s (n = %d)', CHANNELS{c}, cues(k), cueNames{k}, counts(k)));
        xlabel(ax, 'Time [s]');
        ylabel(ax, 'Frequency [Hz]');
    end
end

exportgraphics(fig, OUT_PDF, 'ContentType', 'image', 'Resolution', 300);
fprintf('Figure saved to %s\n', OUT_PDF);

%% Interpretation

% Both hands (773): beta ERD (18-26 Hz) during feedback on C3 and C4,
% stronger on C3 (mean -47%, down to -71%) than on C4 (mean -30%).
% Mu ERD is weaker and Cz barely changes.
% Both feet (771): strongest beta ERD on Cz (mean -51%, down to -84%),
% where the feet area is. C3 shows a large mu ERS (mean +89%).
% Same picture as Lab05, but the maps show the ERD is mostly in 18-26 Hz
% and starts about 1 s after feedback onset.
% Limitations: one subject, 45 trials per class, color scale saturated at
% the 99th percentile, and a few large ERS values pull the means up.

%% Local functions

function lbl = labelEvents(EVENT, sizes)
    n = sum(sizes);
    Tk = zeros(n, 1);
    Fk = zeros(n, 1);
    Ak = zeros(n, 1);
    CFk = zeros(n, 1);
    Rk = zeros(n, 1);

    last = 1;
    ends = cumsum(sizes);
    for i = 1:numel(sizes)
        Rk(last:ends(i)) = i;
        last = ends(i) + 1;
    end

    trial = 0;
    startTrial = NaN;
    for i = 1:numel(EVENT.POS)
        a = EVENT.POS(i);
        b = a + EVENT.DUR(i) - 1;
        if i > 1 && Rk(EVENT.POS(i-1)) ~= Rk(a)
            trial = 0;
        end
        switch EVENT.TYP(i)
            case 1
                trial = trial + 1;
                startTrial = a;
            case 786
                Fk(a:b) = 786;
            case {771, 773}
                Ak(a:b) = EVENT.TYP(i);
            case 781
                CFk(a:b) = 781;
                Tk(startTrial:b) = trial;
            otherwise
                error('Unknown event type %d', EVENT.TYP(i));
        end
    end
    lbl = struct('Tk', Tk, 'Fk', Fk, 'Ak', Ak, 'CFk', CFk, 'Rk', Rk);
end

function [trials, Ck, RunTrial, info] = extractTrials(X, lbl, maxLengthDiff)
    starts = [];
    lens = [];
    Ck = [];
    RunTrial = [];
    TrialInRun = [];
    cueOnset = [];
    feedbackOnset = [];
    fixation = zeros(0, 2);
    feedback = zeros(0, 2);

    for r = unique(lbl.Rk)'
        for t = unique(lbl.Tk(lbl.Rk == r & lbl.Tk > 0))'
            idx = find(lbl.Tk == t & lbl.Rk == r);
            idx = idx(find(lbl.Fk(idx) ~= 0, 1):end);

            cue = unique(lbl.Ak(idx));
            cue(cue == 0) = [];
            assert(numel(cue) == 1, 'Run %d, trial %d: %d cues', r, t, numel(cue));

            fixRel = find(lbl.Fk(idx) ~= 0);
            fbRel = find(lbl.CFk(idx) ~= 0);

            starts(end+1, 1) = idx(1);
            lens(end+1, 1) = numel(idx);
            Ck(end+1, 1) = cue;
            RunTrial(end+1, 1) = r;
            TrialInRun(end+1, 1) = t;
            cueOnset(end+1, 1) = find(lbl.Ak(idx) ~= 0, 1) - 1;
            feedbackOnset(end+1, 1) = fbRel(1) - 1;
            fixation(end+1, :) = [fixRel(1) fixRel(end)];
            feedback(end+1, :) = [fbRel(1) fbRel(end)];
        end
    end

    minLen = min(lens);
    maxLen = max(lens);
    assert(maxLen - minLen <= maxLengthDiff, 'Trial lengths %d..%d windows', minLen, maxLen);
    fixation = min(fixation, minLen);
    feedback = min(feedback, minLen);
    assert(all(fixation(:, 1) <= fixation(:, 2)) && all(feedback(:, 1) <= feedback(:, 2)));

    trials = zeros(minLen, size(X, 2), numel(starts));
    for j = 1:numel(starts)
        trials(:, :, j) = X(starts(j):starts(j) + minLen - 1, :);
    end

    info = struct('length', lens, 'trialInRun', TrialInRun, 'cueOnset', cueOnset, ...
        'feedbackOnset', feedbackOnset, 'fixation', fixation, 'feedback', feedback, ...
        'minLength', minLen, 'maxLength', maxLen);
end
