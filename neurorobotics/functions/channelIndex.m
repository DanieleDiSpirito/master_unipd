function [idx, labels] = channelIndex(names)
% CHANNELINDEX Index of the EEG channels given their 10-20 names.
%
%   [idx, labels] = channelIndex({'C3', 'Cz', 'C4'})   % -> [7 9 11]
%
%   The GDF headers only contain generic labels (eeg:1 ... eeg:16, trigger:1),
%   so the montage is written here. It is consistent with lab05/laplacian16.mat
%   (e.g. Cz neighbours are FCz, C1, C2, CPz). Channel 17 is the trigger.

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
