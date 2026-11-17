%% 1.1

files = { ...
    'ah7.20170613.161402.offline.mi.mi_bhbf.gdf', ...
    'ah7.20170613.162331.offline.mi.mi_bhbf.gdf', ...
    'ah7.20170613.162934.offline.mi.mi_bhbf.gdf'};

s = {0, 0, 0};
h = {0, 0, 0};
sizes = {0, 0, 0};

for i = 1:numel(files)
    [s{i}, h{i}, sizes{i}] = loadRuns(files{i});
end
samplerate = h{1}.SampleRate;
for i = 1:numel(files)
    fprintf('Run %d: %s, %d samples, %d Hz\n', i, files{i}, sizes{i}, h{i}.SampleRate);
end

% sample rate = 512 Hz
% 17 channels
% sizes 172032, 162816, 165888
%% 1.2 - 1.3

% remove trigger, keep the 16 EEG channels in recorded order
for i = 1:numel(files)
    s{i} = s{i}(:, 1:end-1);
end
[~, chLabels] = channelIndex();

L = load('laplacian16.mat');
lap = L.lap;

for i = 1:numel(files)
    assert(size(s{i}, 2) == 16);
end
assert(isequal(size(lap), [16 16]) && all(isfinite(lap(:))));

s_lap = {0, 0, 0};

for i = 1:numel(files)
    s_lap{i} = s{i} * lap;
end

for i = 1:numel(files)
    assert(isequal(size(s_lap{i}), size(s{i})) && all(isfinite(s_lap{i}(:))));
end
%% 2.1

wlength = 0.5; % seconds. Length of the external window
pshift = 0.25; % seconds. Shift of the internal windows
wshift = 0.0625; % seconds. Shift of the external window
mlength = 1; % seconds
PSD = {0, 0, 0};
f = {0, 0, 0};
for i = 1:numel(files)
    [PSD{i}, f{i}] = proc_spectrogram(s_lap{i}, wlength, wshift, pshift, samplerate, mlength);
end

fprintf('Spectrogram: wlength %g s (%d samples), wshift %g s (%d samples), pshift %g s, mlength %g s\n', ...
    wlength, wlength*samplerate, wshift, wshift*samplerate, pshift, mlength);
fprintf('Frequency grid: %g..%g Hz, step %g Hz (%d bins)\n', f{1}(1), f{1}(end), f{1}(2) - f{1}(1), numel(f{1}));

% 2.2
% PSD{1} dim = (windows, frequencies, channels) == (5361, 129, 16)
% PSD{2} dim = (windows, frequencies, channels) == (5073, 129, 16)
% PSD{3} dim = (windows, frequencies, channels) == (5169, 129, 16)
%% 2.3

for i = 1:numel(files)
    PSD{i} = PSD{i}(:, (4:2:48) / 2 + 1, :); % band 4-48Hz selected
    f{i} = f{i}((4:2:48) / 2 + 1);
end
fprintf('Selected frequencies: %g..%g Hz, step 2 Hz (%d bins)\n', f{1}(1), f{1}(end), numel(f{1}));
%% 2.check

for i = 1:numel(files)
    assert(all(isfinite(PSD{i}(:))) && all(PSD{i}(:) >= 0));
    assert(numel(f{i}) == size(PSD{i}, 2));
    assert(size(PSD{i}, 3) == numel(chLabels));
end
%% 3.1 - 3.2

events = {0, 0, 0};
for i = 1:numel(files)
    events{i} = h{i}.EVENT;
end

POSwin = {0, 0, 0};
ENDwin = {0, 0, 0};
DURwin = {0, 0, 0};
winconv = 'backward';
for i = 1:numel(files)
    POSwin{i} = proc_pos2win(events{i}.POS, wshift*samplerate, winconv, wlength*samplerate);
    ENDwin{i} = proc_pos2win(events{i}.POS + events{i}.DUR - 1, wshift*samplerate, winconv, wlength*samplerate);
    DURwin{i} = ENDwin{i} - POSwin{i} + 1;
end
fprintf('Event conversion: proc_pos2win, ''%s'', wPOS = floor(POS/%d) + 1 - %d; end = POS + DUR - 1\n', ...
    winconv, wshift*samplerate, floor(wlength*samplerate / (wshift*samplerate)));
%% 3.3

% Events from proc_pos2win are ~0.44 s (7 windows) late: the helper ignores
% the 1 s moving average of proc_spectrogram (checked with a test signal).
for i = 1:numel(files)
    assert(sum(POSwin{i} < 1) == 0);
    assert(sum(ENDwin{i} > size(PSD{i}, 1)) == 0);
    assert(sum(DURwin{i} <= 0) == 0);
    fprintf('Run %d: first event sample %d -> window %d, last event ends at window %d of %d\n', ...
        i, events{i}.POS(1), POSwin{i}(1), ENDwin{i}(end), size(PSD{i}, 1));
end
%% 3.4

for i = 1:numel(files)
    [~, base] = fileparts(files{i});
    out = struct();
    out.PSD = PSD{i};
    out.f = f{i};
    out.chLabels = chLabels;
    out.samplerate = samplerate;
    out.wlength = wlength;
    out.wshift = wshift;
    out.pshift = pshift;
    out.mlength = mlength;
    out.winconv = winconv;
    out.EVENT = events{i};
    out.EVENTwin = struct('POS', POSwin{i}, 'DUR', DURwin{i}, 'TYP', events{i}.TYP);
    out.run = i;
    out.nSamples = size(s{i}, 1);
    out.nWindows = size(PSD{i}, 1);
    save([base '.mat'], '-struct', 'out');
    fprintf('Saved %s.mat (run %d, %d windows)\n', base, i, out.nWindows);
end

%% Local functions

function [s, h, sizes] = loadRuns(files)
% Load the GDF runs, check sample rate and channels, concatenate and shift event positions.
    if ischar(files) || isstring(files)
        files = cellstr(files);
    end
    N_RUNS = numel(files);
    ss = cell(N_RUNS, 1);
    hs = cell(N_RUNS, 1);
    for i = 1:N_RUNS
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
