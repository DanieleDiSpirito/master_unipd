function plotEvents(ax, h, fs, tLim, yLim)
% PLOTEVENTS Draw the events falling in [tLim(1), tLim(2)] seconds as colored patches.
%
%   plotEvents(ax, h, fs, [START END], [-80 80])
%
%   Time of sample n is n/fs (same convention used to plot the signal).
%   Adds a legend with only the event codes visible in the window.

    colors = eventColors();
    hold(ax, 'on');

    shown = [];
    for j = 1:numel(h.EVENT.POS)
        start = h.EVENT.POS(j) / fs;
        term  = start + h.EVENT.DUR(j) / fs;
        code  = h.EVENT.TYP(j);

        if term < tLim(1) || start > tLim(2)
            continue
        end

        patch(ax, [start term term start], [yLim(1) yLim(1) yLim(2) yLim(2)], ...
            colors(code), ...
            'FaceAlpha', 0.25, ...
            'EdgeColor', 'none', ...
            'HandleVisibility', 'off');

        xline(ax, start, 'k-', 'HandleVisibility', 'off');
        xline(ax, term, 'k-', 'HandleVisibility', 'off');
        shown(end+1) = code; %#ok<AGROW>
    end

    % One dummy patch per visible code, only for the legend
    codes = unique(shown);
    p = gobjects(numel(codes), 1);
    for k = 1:numel(codes)
        p(k) = patch(ax, NaN, NaN, colors(codes(k)), ...
            'FaceAlpha', 0.25, 'EdgeColor', 'none', ...
            'DisplayName', int2str(codes(k)));
    end
    if ~isempty(p)
        legend(p, 'Location', 'bestoutside');
    end

    xlim(ax, tLim);
    ylim(ax, yLim);
end
