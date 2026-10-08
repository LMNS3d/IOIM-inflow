function ha = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
%TIGHT_SUBPLOT Create subplot axes with explicit gaps and margins.
%   HA = TIGHT_SUBPLOT(NH, NW, GAP, MARG_H, MARG_W) fills the current
%   figure with NH-by-NW axes and returns their handles in row-major
%   order.
%
%   GAP    [gap_vertical, gap_horizontal], normalized figure units.
%   MARG_H [bottom, top].
%   MARG_W [left, right].
%
%   Call signature matches Pekka Kumpulainen's File Exchange utility
%   (id 27991). This file is an independent reimplementation so the
%   repository does not redistribute that third-party source. Use it
%   only for figure layout; it does not affect the flow solution.
%
%   See also SUBPLOT, TILEDLAYOUT.

if nargin < 3 || isempty(gap),    gap = [0.02, 0.02]; end
if nargin < 4 || isempty(marg_h), marg_h = [0.05, 0.05]; end
if nargin < 5 || isempty(marg_w), marg_w = [0.05, 0.05]; end

if numel(gap) == 1,    gap = [gap, gap]; end
if numel(marg_h) == 1, marg_h = [marg_h, marg_h]; end
if numel(marg_w) == 1, marg_w = [marg_w, marg_w]; end

axh = (1 - sum(marg_h) - (Nh - 1) * gap(1)) / Nh;
axw = (1 - sum(marg_w) - (Nw - 1) * gap(2)) / Nw;

if axh <= 0 || axw <= 0
    error('tight_subplot:BadMargins', ...
        'Gaps and margins leave no room for the axes.');
end

ha = gobjects(Nh * Nw, 1);
ii = 0;
for ih = 1:Nh
    for iw = 1:Nw
        ii = ii + 1;
        x = marg_w(1) + (iw - 1) * (axw + gap(2));
        y = 1 - marg_h(2) - ih * axh - (ih - 1) * gap(1);
        ha(ii) = axes('Units', 'normalized', 'Position', [x, y, axw, axh]);
    end
end
end
