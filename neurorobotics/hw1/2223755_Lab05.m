clear; close all; clc;

%% Configuration
GDF_FILES = { ...
    'ah7.20170613.161402.offline.mi.mi_bhbf.gdf', ...
    'ah7.20170613.162331.offline.mi.mi_bhbf.gdf', ...
    'ah7.20170613.162934.offline.mi.mi_bhbf.gdf'};
LAP_FILE = 'laplacian16.mat';
CHANLOCS_FILE = 'chanlocs16.mat';
BIOSIG_DIR = '../toolboxes/biosig/biosig';
EEGLAB_DIR = '../toolboxes/eeglab';
OUT_PDF = 'StudentID_Lab05.pdf';

cues = [773 771]; % hands, feet (order of the topoplots)
cueNames = {'both hands', 'both feet'};
cueColors = [0.20 0.70 0.35; 0.20 0.60 0.90];
CHANNEL = 'C3';

muButter = [8 12];   % Hz
betaButter = [18 22]; % Hz
ORDER = 3;            % Butterworth order (band-pass -> 2*ORDER poles), filtfilt
WINDOW = 1;           % s, moving average
MAX_LENGTH_DIFF = 5;  % samples, event-duration jitter between trials (<= 10 ms)

%% Toolboxes
if isempty(which('sload'))
    p = strsplit(genpath(BIOSIG_DIR), pathsep);
    addpath(strjoin(p(~cellfun(@isempty, p) & ~contains(p, 'maybe-missing')), pathsep));
end
if isempty(which('topoplot'))
    addpath(EEGLAB_DIR);
    eeglab nogui;
end

%% 1. Loading, events, Laplacian
[s, h, sizes] = loadRuns(GDF_FILES);
SampleRate = h.SampleRate;

lbl = labelEvents(h, sizes);
assert(all(structfun(@numel, lbl) == size(s, 1)));
assert(isequal(unique(lbl.Ak(lbl.Ak > 0)), sort(cues(:))));

s = s(:, 1:end-1); % remove trigger, 16 EEG channels in recorded order
[~, chLabels] = channelIndex();

L = load(LAP_FILE);
lap = L.lap;
assert(size(s, 2) == 16);
assert(isequal(size(lap), [16 16]) && all(isfinite(lap(:))));

s_lap = s * lap;
assert(isequal(size(s_lap), size(s)) && all(isfinite(s_lap(:))));

%% 2. Band-pass filtering and smoothed power (no log), per run
sfiltmu = bandpassRuns(s_lap, lbl.Rk, muButter, ORDER, SampleRate);
sfiltbe = bandpassRuns(s_lap, lbl.Rk, betaButter, ORDER, SampleRate);

P_mu = bandPower(sfiltmu, lbl.Rk, SampleRate, WINDOW);
P_be = bandPower(sfiltbe, lbl.Rk, SampleRate, WINDOW);

fprintf('Sampling frequency: %d Hz\n', SampleRate);
fprintf('Butterworth order %d, filtfilt: mu [%g %g] Hz, beta [%g %g] Hz\n', ORDER, muButter, betaButter);
fprintf('Moving average: %g s = %d samples (filter, per run)\n', WINDOW, round(WINDOW * SampleRate));

%% 3. Trials (fixation onset -> end of feedback) and periods
[trials, Ck, RunTrial, info] = extractTrials({P_mu, P_be}, lbl, MAX_LENGTH_DIFF);
[TrialPower_mu, TrialPower_be] = trials{:};
assert(size(unique([RunTrial info.trialInRun], 'rows'), 1) == numel(Ck));

refPeriod = info.fixation; % [trials x 2] indices inside each trial
actPeriod = info.feedback;

fprintf('Trials: %d, length %d..%d -> %d samples\n', numel(Ck), info.minLength, info.maxLength, size(TrialPower_mu, 1));
fprintf('Reference = fixation, samples %d..%d (end %d..%d)\n', min(refPeriod(:, 1)), max(refPeriod(:, 1)), min(refPeriod(:, 2)), max(refPeriod(:, 2)));
fprintf('Activity = continuous feedback, samples %d..%d (end %d..%d)\n', min(actPeriod(:, 1)), max(actPeriod(:, 1)), min(actPeriod(:, 2)), max(actPeriod(:, 2)));

%% 4. ERD/ERS = 100 * (P - R) ./ R, R = fixation mean per channel and trial
[ERD_ERS_mu, ~, valid_mu] = erdErs(TrialPower_mu, refPeriod);
[ERD_ERS_be, ~, valid_be] = erdErs(TrialPower_be, refPeriod);

valid = valid_mu & valid_be;
fprintf('Excluded trials: %d\n', nnz(~valid));
for j = find(~valid)'
    fprintf(' - run %d, trial %d (cue %d)\n', RunTrial(j), info.trialInRun(j), Ck(j));
end

ERD_ERS = {ERD_ERS_mu(:, :, valid), ERD_ERS_be(:, :, valid)};
Ck = Ck(valid);
refPeriod = refPeriod(valid, :);
actPeriod = actPeriod(valid, :);
cueOn = info.cueOnset(valid);
fbOn = info.feedbackOnset(valid);

for b = 1:numel(ERD_ERS)
    assert(all(abs(periodMean(ERD_ERS{b}, refPeriod)) < 1e-6, 'all')); % fixation mean ~0
end

fprintf('Cue onset %.3f..%.3f s, feedback %.3f..%.3f s -> %.3f s (from fixation)\n', ...
    min(cueOn) / SampleRate, max(cueOn) / SampleRate, min(fbOn) / SampleRate, max(fbOn) / SampleRate, ...
    max(actPeriod(:, 2) - 1) / SampleRate);

%% 5. Temporal averages and activity maps
channel = channelIndex(CHANNEL);
bandNames = {sprintf('\\mu [%g-%g] Hz', muButter), sprintf('\\beta [%g-%g] Hz', betaButter)};

avg = cell(numel(ERD_ERS), 1);
sem = cell(numel(ERD_ERS), 1);
actMap = cell(numel(ERD_ERS), 1); % {band} -> [classes x channels]
for b = 1:numel(ERD_ERS)
    [avg{b}, counts, ~, sem{b}] = cueAverage(ERD_ERS{b}(:, channel, :), Ck, cues);
    act = periodMean(ERD_ERS{b}, actPeriod);
    for c = 1:numel(cues)
        actMap{b}(c, :) = mean(act(Ck == cues(c), :), 1);
    end
end

C = load(CHANLOCS_FILE);
chanlocs = matchChanlocs(C.chanlocs16, chLabels);

%% 6. Figure
t = (0:size(ERD_ERS{1}, 1)-1) ./ SampleRate; % time from fixation onset
cmap = blueWhiteRed();

fig = figure('Color', 'w', 'Position', [50 50 1600 950]);

for b = 1:numel(ERD_ERS)
    ax = subplot(2, 4, (b - 1) * 2 + [1 2]);
    hold(ax, 'on');
    for c = 1:numel(cues)
        m = avg{b}(:, 1, c);
        e = sem{b}(:, 1, c);
        fill(ax, [t fliplr(t)], [m + e; flipud(m - e)]', cueColors(c, :), ...
            'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        plot(ax, t, m, 'Color', cueColors(c, :), 'LineWidth', 1.5, ...
            'DisplayName', sprintf('%d %s (n = %d)', cues(c), cueNames{c}, counts(c)));
    end
    yline(ax, 0, 'k-', 'HandleVisibility', 'off');
    markTrialEvents(ax, cueOn, fbOn, SampleRate, actPeriod(:, 2) - 1);
    hold(ax, 'off');
    title(ax, sprintf('%s - %s (mean \\pm 1 SEM)', bandNames{b}, chLabels{channel}));
    xlabel(ax, 'time from fixation onset [s]');
    ylabel(ax, 'ERD/ERS [%]');
    legend(ax, 'show', 'Location', 'southwest');
    xlim(ax, [0, t(end)]);
    grid(ax, 'on');
end

for b = 1:numel(ERD_ERS)
    lim = max(abs(actMap{b}(:))); % same symmetric limits for both classes of a band
    for c = 1:numel(cues)
        ax = subplot(2, 4, 4 + (b - 1) * numel(cues) + c);
        topoplot(actMap{b}(c, :), chanlocs, 'maplimits', [-lim lim], 'colormap', cmap, 'electrodes', 'labels');
        colormap(ax, cmap);
        clim(ax, [-lim lim]);
        cb = colorbar(ax, 'southoutside');
        cb.Label.String = 'ERD/ERS (%)';
        title(ax, sprintf('%s - %s\ncontinuous feedback (n = %d)', bandNames{b}, cueNames{c}, counts(c)));
    end
end

sgtitle(fig, sprintf(['Lab05 - Laplacian EEG, 3 offline runs - reference: fixation, activity: continuous feedback\n' ...
    'Butterworth order %d + filtfilt, %g s moving average - 773 = both hands, 771 = both feet'], ORDER, WINDOW));

% topoplot (EEGLAB) uses a light-blue background for the figure and the head ring
set(fig, 'Color', 'w');
set(findobj(fig, 'Type', 'patch', 'FaceColor', [.93 .96 1]), 'FaceColor', 'w');
exportgraphics(fig, OUT_PDF, 'ContentType', 'image', 'Resolution', 300, 'BackgroundColor', 'white');
fprintf('Figure saved to %s\n', OUT_PDF);

%% Interpretation
% C3, feedback period (45 trials per class):
% - mu: opposite sign between classes. Hands -> ERD (mean -16%, median -34%),
%   feet -> ERS (mean +32%, up to ~+60% at the end of feedback).
% - beta: ERD for both, stronger and more consistent for hands (mean -53%,
%   min ~-68% at ~5.8 s) than for feet (mean -20%, median -56%).
% - The separation builds up ~1 s after feedback onset (partly the delay of
%   the 1 s moving average); around the cue the pattern differs.
% Maps:
% - beta hands: bilateral ERD, strongest on C3 (-53%) and C2/C4 (~-38%).
% - beta feet: ERD strongest on Cz (-41%), but also on Fz/FC3/CP3 (less focal).
% - mu feet: ERD on CPz/Cz with ERS on C1/C3 (focal ERD, surround ERS).
% - mu hands: weak ERD on C3/C1 and a strong ERS on CPz (+99%, median +67%),
%   not the expected pattern.
% Limitations: one subject, 45 trials per class; SEM bands are not a test;
% ERD/ERS is skewed (a few trials up to +622% move the means); effects vary
% across runs; event timing varies <= 10 ms; map colors far from the 16
% central electrodes are interpolated.

%% Local functions

function [s, h, sizes] = loadRuns(files)
% Load the GDF runs, check sample rate and channels, concatenate and shift event positions.
    N_RUNS = numel(files);
    ss = cell(N_RUNS, 1);
    hs = cell(N_RUNS, 1);
    for i = 1:N_RUNS
        assert(isfile(files{i}), 'File not found: %s', files{i});
        [ss{i}, hs{i}] = sload(files{i});
    end

    h = hs{1};
    for i = 2:N_RUNS
        assert(hs{i}.SampleRate == h.SampleRate, 'Run %d: different sample rate', i);
        assert(isequal(cellstr(hs{i}.Label), cellstr(h.Label)), 'Run %d: different channels', i);
    end

    sizes = cellfun(@(x) size(x, 1), ss);
    offsets = [0; cumsum(sizes(1:end-1))];
    s = cat(1, ss{:});

    POS = [];
    TYP = [];
    DUR = [];
    for i = 1:N_RUNS
        POS = [POS; hs{i}.EVENT.POS + offsets(i)]; %#ok<AGROW>
        TYP = [TYP; hs{i}.EVENT.TYP]; %#ok<AGROW>
        DUR = [DUR; hs{i}.EVENT.DUR]; %#ok<AGROW>
    end
    h.EVENT = struct('POS', POS, 'TYP', TYP, 'DUR', DUR);

    assert(size(s, 1) == sum(sizes));
    assert(numel(POS) == sum(cellfun(@(x) numel(x.EVENT.POS), hs)));
end

function lbl = labelEvents(h, sizes)
% Sample-wise labels: Tk (trial inside the run), Fk, Ak, CFk, Xk, Rk.
    nSamples = sum(sizes);
    Tk = zeros(nSamples, 1);
    Fk = zeros(nSamples, 1);
    Ak = zeros(nSamples, 1);
    CFk = zeros(nSamples, 1);
    Xk = zeros(nSamples, 1);
    Rk = zeros(nSamples, 1);

    last_pos = 1;
    offsets = cumsum(sizes);
    for i = 1:numel(sizes)
        Rk(last_pos:offsets(i)) = i;
        last_pos = offsets(i) + 1;
    end

    POS = h.EVENT.POS;
    TYP = h.EVENT.TYP;
    DUR = h.EVENT.DUR;

    trial = 0;
    trialOpen = false;
    startTrialIdx = NaN;
    for i = 1:numel(POS)
        startIdx = POS(i);
        endIdx = startIdx + DUR(i) - 1;
        assert(Rk(startIdx) == Rk(max(endIdx, startIdx)), 'Event %d crosses a run boundary', i);

        if i > 1 && Rk(POS(i-1)) ~= Rk(startIdx)
            trial = 0;
            trialOpen = false;
        end

        switch TYP(i)
            case 1 % start
                trial = trial + 1;
                startTrialIdx = startIdx;
                trialOpen = true;
            case 786 % fixation
                if ~trialOpen
                    trial = trial + 1;
                    startTrialIdx = startIdx;
                    trialOpen = true;
                end
                Fk(startIdx:endIdx) = TYP(i);
            case {771, 773, 783} % cue
                Ak(startIdx:endIdx) = TYP(i);
            case 781 % continuous feedback
                CFk(startIdx:endIdx) = TYP(i);
                Tk(startTrialIdx:endIdx) = trial;
                trialOpen = false;
            case {897, 898} % hit / miss
                Xk(startIdx:endIdx) = TYP(i);
            otherwise
                error('Weird stuff in GDFs: event type %d', TYP(i));
        end
    end

    lbl = struct('Tk', Tk, 'Fk', Fk, 'Ak', Ak, 'CFk', CFk, 'Xk', Xk, 'Rk', Rk);
end

function [trials, Ck, RunTrial, info] = extractTrials(X, lbl, maxLengthDiff)
% Trials from fixation onset to end of feedback, selected with BOTH Tk and Rk,
% trimmed at the END to the shortest length. X is a cell of [samples x channels].
    starts = [];
    lens = [];
    Ck = [];
    RunTrial = [];
    TrialInRun = [];
    cueOnset = [];
    feedbackOnset = [];
    fixation = zeros(0, 2);
    feedback = zeros(0, 2);

    for r = unique(lbl.Rk(lbl.Rk > 0))'
        for t = unique(lbl.Tk(lbl.Rk == r & lbl.Tk > 0))'
            idx = find(lbl.Tk == t & lbl.Rk == r);
            assert(all(diff(idx) == 1), 'Trial %d of run %d is not contiguous', t, r);
            idx = idx(find(lbl.Fk(idx) ~= 0, 1):end); % start at fixation

            cue = unique(lbl.Ak(idx));
            cue(cue == 0) = [];
            assert(numel(cue) == 1, 'Trial %d of run %d: %d cues', t, r, numel(cue));

            fixRel = find(lbl.Fk(idx) ~= 0);
            fbRel = find(lbl.CFk(idx) ~= 0);
            assert(~isempty(fixRel) && ~isempty(fbRel), 'Trial %d of run %d: empty period', t, r);

            starts(end+1, 1) = idx(1); %#ok<AGROW>
            lens(end+1, 1) = numel(idx); %#ok<AGROW>
            Ck(end+1, 1) = cue; %#ok<AGROW>
            RunTrial(end+1, 1) = r; %#ok<AGROW>
            TrialInRun(end+1, 1) = t; %#ok<AGROW>
            cueOnset(end+1, 1) = find(lbl.Ak(idx) ~= 0, 1) - 1; %#ok<AGROW>
            feedbackOnset(end+1, 1) = fbRel(1) - 1; %#ok<AGROW>
            fixation(end+1, :) = [fixRel(1) fixRel(end)]; %#ok<AGROW>
            feedback(end+1, :) = [fbRel(1) fbRel(end)]; %#ok<AGROW>
        end
    end

    minLen = min(lens);
    maxLen = max(lens);
    assert(maxLen - minLen <= maxLengthDiff, 'Trial lengths %d..%d: check the extraction', minLen, maxLen);

    fixation = min(fixation, minLen);
    feedback = min(feedback, minLen);
    assert(all(fixation(:, 1) <= fixation(:, 2)) && all(feedback(:, 1) <= feedback(:, 2)));

    trials = cell(size(X));
    for k = 1:numel(X)
        trials{k} = zeros(minLen, size(X{k}, 2), numel(starts));
        for j = 1:numel(starts)
            trials{k}(:, :, j) = X{k}(starts(j):starts(j) + minLen - 1, :);
        end
    end

    info = struct('length', lens, 'trialInRun', TrialInRun, 'cueOnset', cueOnset, ...
        'feedbackOnset', feedbackOnset, 'fixation', fixation, 'feedback', feedback, ...
        'minLength', minLen, 'maxLength', maxLen);
end

function sfilt = bandpassRuns(s, Rk, band, order, fs)
% Zero-phase Butterworth band-pass, each run filtered separately.
    assert(band(1) > 0 && band(1) < band(2) && band(2) < fs / 2, 'Invalid cutoffs');
    [b, a] = butter(order, 2 * band / fs);
    sfilt = zeros(size(s));
    for r = unique(Rk(Rk > 0))'
        idx = Rk == r;
        sfilt(idx, :) = filtfilt(b, a, s(idx, :));
    end
    assert(all(isfinite(sfilt(:))));
end

function P = bandPower(sfilt, Rk, fs, winSec)
% Square, then trailing moving average with filter (state reset at each run).
    windowSize = round(winSec * fs);
    b = ones(1, windowSize) / windowSize;
    a = 1;
    sq = sfilt .^ 2;
    P = zeros(size(sfilt));
    for r = unique(Rk(Rk > 0))'
        idx = Rk == r;
        P(idx, :) = filter(b, a, sq(idx, :), [], 1);
    end
    assert(all(isfinite(P(:))) && all(P(:) >= 0), 'Power must be finite and non-negative');
end

function [ERD_ERS, Reference, valid] = erdErs(TrialPower, refPeriod)
% Percentage change w.r.t. the mean power of each trial's own reference period.
    nTrials = size(TrialPower, 3);
    Reference = zeros(1, size(TrialPower, 2), nTrials);
    for j = 1:nTrials
        Reference(1, :, j) = mean(TrialPower(refPeriod(j, 1):refPeriod(j, 2), :, j), 1);
    end
    valid = squeeze(all(isfinite(Reference) & Reference > 0, 2));
    valid = valid(:);
    ERD_ERS = 100 * (TrialPower - Reference) ./ Reference;
    ERD_ERS(:, :, ~valid) = NaN;
end

function [avg, counts, cues, sem] = cueAverage(trials, Ck, cues)
% Mean and standard error (sample std / sqrt(n)) across the trials of each cue.
    avg = zeros(size(trials, 1), size(trials, 2), numel(cues));
    sem = zeros(size(avg));
    counts = zeros(numel(cues), 1);
    for i = 1:numel(cues)
        sel = Ck == cues(i);
        counts(i) = nnz(sel);
        avg(:, :, i) = mean(trials(:, :, sel), 3);
        sem(:, :, i) = std(trials(:, :, sel), 0, 3) / sqrt(counts(i));
    end
end

function values = periodMean(trials, period)
% [trials x channels] mean of each trial over its own period.
    values = zeros(size(trials, 3), size(trials, 2));
    for j = 1:size(trials, 3)
        values(j, :) = mean(trials(period(j, 1):period(j, 2), :, j), 1);
    end
end

function [idx, labels] = channelIndex(names)
% 16-channel montage (GDF headers only contain eeg:1..16).
    labels = {'Fz', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', ...
              'C3', 'C1', 'Cz', 'C2', 'C4', ...
              'CP3', 'CP1', 'CPz', 'CP2', 'CP4'};
    if nargin < 1
        idx = 1:numel(labels);
        return
    end
    names = cellstr(names);
    [found, idx] = ismember(lower(names), lower(labels));
    assert(all(found), 'Unknown channel(s): %s', strjoin(names(~found), ', '));
end

function chanlocs = matchChanlocs(chanlocs, chLabels)
% Reorder the locations by label to match the EEG channels.
    locLabels = strtrim({chanlocs.labels});
    assert(numel(unique(lower(locLabels))) == numel(locLabels), 'Duplicated labels in chanlocs');
    [found, order] = ismember(lower(chLabels), lower(locLabels));
    assert(all(found), 'Missing location for: %s', strjoin(chLabels(~found), ', '));
    chanlocs = chanlocs(order);
    assert(isequal(lower(strtrim({chanlocs.labels})), lower(chLabels(:)')));
end

function markTrialEvents(ax, cueOnset, feedbackOnset, fs, feedbackEnd)
% Median event lines; shaded [min, max] range when onsets differ across trials.
    names = {'Cue', 'Feedback', 'Feedback end'};
    onsets = {cueOnset(:) / fs, feedbackOnset(:) / fs, feedbackEnd(:) / fs};
    hold(ax, 'on');
    for k = 1:numel(onsets)
        o = onsets{k};
        if max(o) > min(o)
            xregion(ax, min(o), max(o), 'FaceColor', [0.5 0.5 0.5], 'HandleVisibility', 'off');
        end
        xline(ax, median(o), 'k--', names{k}, 'HandleVisibility', 'off');
    end
end

function cmap = blueWhiteRed(n)
% Diverging colormap blue - white - red (use symmetric color limits).
    if nargin < 1
        n = 256;
    end
    half = floor(n / 2);
    x = linspace(0, 1, half)';
    cmap = [x x ones(half, 1); ones(half, 1) flipud(x) flipud(x)];
end
