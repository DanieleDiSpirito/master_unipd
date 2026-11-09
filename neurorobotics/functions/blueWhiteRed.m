function cmap = blueWhiteRed(n)
% BLUEWHITERED Diverging colormap blue (negative) - white (zero) - red (positive).
%
%   cmap = blueWhiteRed(n)   % n colors, default 256
%
%   Centered at zero only if the color limits are symmetric, e.g. clim([-L L]).

    if nargin < 1
        n = 256;
    end
    half = floor(n / 2);
    x = linspace(0, 1, half)';
    blue = [x x ones(half, 1)];               % blue -> white
    red = [ones(half, 1) flipud(x) flipud(x)]; % white -> red
    cmap = [blue; red];
    if mod(n, 2) == 1
        cmap = [blue; 1 1 1; red];
    end
end
