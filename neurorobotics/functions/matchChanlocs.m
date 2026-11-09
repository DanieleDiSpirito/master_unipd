function chanlocs = matchChanlocs(chanlocs, chLabels)
% MATCHCHANLOCS Reorder channel locations to match the EEG channel labels.
%
%   chanlocs = matchChanlocs(chanlocs, chLabels)
%
%   chanlocs : EEGLAB chanlocs struct array (field .labels)
%   chLabels : {1 x channels} EEG labels in recorded order (e.g. channelIndex)
%
%   Output chanlocs(j) is the location of chLabels{j}. Errors if a label is
%   missing or duplicated (same dimensions do NOT imply same order).

    locLabels = strtrim({chanlocs.labels});
    assert(numel(unique(lower(locLabels))) == numel(locLabels), 'Duplicated labels in chanlocs');

    [found, order] = ismember(lower(chLabels), lower(locLabels));
    assert(all(found), 'Missing location for: %s', strjoin(chLabels(~found), ', '));

    if ~isequal(order(:)', 1:numel(chLabels))
        fprintf('matchChanlocs: chanlocs reordered\n');
    end
    chanlocs = chanlocs(order);

    assert(isequal(lower(strtrim({chanlocs.labels})), lower(chLabels(:)')));
end
