function [s, h, sizes, hs] = loadRuns(files)
% LOADRUNS Load one or more GDF runs, check they are compatible and concatenate them.
%
%   [s, h, sizes, hs] = loadRuns(files)
%
%   files : path (char/string) or cell array of paths, in concatenation order
%
%   s     : [samples x channels] concatenated EEG (all channels, trigger included)
%   h     : header of the first run, with h.EVENT.POS/TYP/DUR of ALL runs
%           (POS shifted by the number of samples already appended)
%   sizes : [runs x 1] number of samples of each run
%   hs    : {runs x 1} original headers (useful for validation)

    if ischar(files) || isstring(files)
        files = cellstr(files);
    end

    N_RUNS = numel(files);
    ss = cell(N_RUNS, 1);
    hs = cell(N_RUNS, 1);

    for i = 1:N_RUNS
        [ss{i}, hs{i}] = sload(files{i});
    end

    % All runs must have the same sampling frequency and channels (same order)
    h = hs{1};
    for i = 2:N_RUNS
        assert(hs{i}.SampleRate == h.SampleRate, ...
            'Run %d has sample rate %d Hz, run 1 has %d Hz', i, hs{i}.SampleRate, h.SampleRate);
        assert(isequal(cellstr(hs{i}.Label), cellstr(h.Label)), ...
            'Run %d has different channels (or order) than run 1', i);
    end

    sizes = cellfun(@(x) size(x, 1), ss);
    offsets = [0; cumsum(sizes(1:end-1))]; % samples already appended before each run

    s = cat(1, ss{:});

    POS = [];
    TYP = [];
    DUR = [];
    for i = 1:N_RUNS
        POS = [POS; hs{i}.EVENT.POS + offsets(i)]; %#ok<AGROW> % no extra sample added
        TYP = [TYP; hs{i}.EVENT.TYP]; %#ok<AGROW>
        DUR = [DUR; hs{i}.EVENT.DUR]; %#ok<AGROW>
    end

    h.EVENT.POS = POS;
    h.EVENT.TYP = TYP;
    h.EVENT.DUR = DUR;
    if isfield(h.EVENT, 'CHN')
        h.EVENT = rmfield(h.EVENT, 'CHN'); % would refer only to run 1
    end

    % Combined counts must equal the sums of the original ones
    assert(size(s, 1) == sum(sizes));
    assert(numel(h.EVENT.POS) == sum(cellfun(@(x) numel(x.EVENT.POS), hs)));
end
