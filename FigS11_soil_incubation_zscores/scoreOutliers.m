% Input: x and y that are expected to hug 1:1 line
% !! Code assumes these are biomass values with pseudocount 0.5 !!
% Assume error model: a constant fractional noise + a constant #counts
% assign z-score to each point (array same size as x)
% Use the params supplied; if none, estimate them from x & y
function [score, cFrac, c0] = scoreOutliers(x, y, cFrac, c0) 
sz = size(x);

% x = lambda xRaw
% min(x) = lambda min(xRaw) = lambda*0.5
% xRaw = x/lambda;
xRaw = x./(min(x)/0.5); 
yRaw = y./(min(y)/0.5); 
x = x(:); y = y(:); xRaw = xRaw(:); yRaw = yRaw(:);

if nargin<3
    % Estimate the fractional noise by looking at the high-abundance ASVs only
    RAW_CUTOFF = 50;
    sel = xRaw>RAW_CUTOFF & yRaw>RAW_CUTOFF;
    ratios = x(sel)./y(sel);
    % For a gaussian:
    % sqrt(mean((ratios-1).^2)) = std
    % sqrt(median((ratios-1).^2)) = 0.67*std
    cFrac = sqrt(median((ratios-1).^2))/0.67;
end

s = 1./mean([min(x)/0.5, min(y)/0.5]); % roughly, the conversion factor from raw counts to values in x and y

% for each point, compute the z-score (actual deviation divided by sigma), for a given c0
% Use "max" because of the two quantities, one is positive, and the other (negative) is irrelevant
howFarOut = @(c0) sign(y-x).*max((y-x)./sqrt((cFrac*x).^2+(c0/s)^2), (x-y)./sqrt((cFrac*y).^2+(c0/s)^2));

if nargin<4
    % Now, keeping cFrac as fixed, estimate c0 as the value for whcih 67%
    % of points are within 1 sigma by this error model

    % Only pay attention to taxa that are indeed present.
    % 1 is safe to use even with rounding, because with pseudocount 0.5,
    % the lowest value for taxa that were present is actually 1.5 
    restrictToPresentTaxaOnly = @(v) v(xRaw(:)>=1 | yRaw(:)>=1);
    fracInside = @(c0)mean( restrictToPresentTaxaOnly(abs(howFarOut(c0))<1) );
    % adjust c0 so the fraction inside is 67% as one expects for a "1-sigma" deviation    
    c0 = fzero(@(c0)fracInside(c0)-0.67, 10);
end
score = howFarOut(c0);
score = reshape(score,sz);
end