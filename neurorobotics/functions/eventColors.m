function colors = eventColors()
% EVENTCOLORS Map event code -> RGB color used in the plots.

    colors = containers.Map('KeyType', 'double', 'ValueType', 'any');

    colors(001) = [0.90 0.20 0.20]; % red        trial start
    colors(771) = [0.20 0.60 0.90]; % sky blue   cue both feet
    colors(773) = [0.20 0.70 0.35]; % green      cue both hands
    colors(781) = [0.65 0.35 0.80]; % purple     continuous feedback
    colors(783) = [0.95 0.60 0.15]; % orange     cue rest
    colors(786) = [0.15 0.75 0.75]; % turquoise  fixation
    colors(897) = [0.80 0.45 0.55]; % pink       hit
    colors(898) = [0.55 0.45 0.25]; % brown      miss
end
