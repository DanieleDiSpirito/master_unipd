function markTrialEvents(ax, cueOnset, feedbackOnset, fs, label, feedbackEnd)
% MARKTRIALEVENTS Mark cue and feedback onsets (samples from trial start) on a trial plot.
%
%   markTrialEvents(ax, info.cueOnset(j), info.feedbackOnset(j), fs, label)   % one trial
%   markTrialEvents(ax, info.cueOnset(sel), info.feedbackOnset(sel), fs, label) % many trials
%   markTrialEvents(..., label, feedbackEnd)  % also mark the end of feedback
%
%   With many trials a line is drawn at the median onset and, if the onsets
%   are not identical, the [min, max] range is shaded.

    names = {'Cue', 'Feedback'};
    onsets = {cueOnset(:) / fs, feedbackOnset(:) / fs};
    if nargin >= 6
        names{end+1} = 'Feedback end';
        onsets{end+1} = feedbackEnd(:) / fs;
    end

    hold(ax, 'on');
    for k = 1:numel(onsets)
        o = onsets{k};
        if max(o) > min(o)
            xregion(ax, min(o), max(o), 'FaceColor', [0.5 0.5 0.5], 'HandleVisibility', 'off');
        end
        if label
            xline(ax, median(o), 'k--', names{k}, 'HandleVisibility', 'off');
        else
            xline(ax, median(o), 'k--', 'HandleVisibility', 'off');
        end
    end
end
