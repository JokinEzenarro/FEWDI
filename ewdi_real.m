% =========================================================================
% ewdi_real.m
%
% Author: Jokin Ezenarro
% Contact: jokin@food.ku.dk
% Last modified: 14/4/2026
%
% License/usage:
% Free to use and modify for non-commercial purposes, provided that proper
% citation is given to the original work and repository.
% =========================================================================

function EWDIr = ewdi_real(sample, loading0, t0)
% EWDI_REAL Compute the measured/real EWDI trajectory for one sample.
%
% INPUTS
%   sample   : (nTime x nVariables) spectral trajectory of one fermentation
%   loading0 : reference loading vector obtained with EWDI_REFERENCE
%   t0       : minimum number of initial points before EWDI is calculated
%
% OUTPUT
%   EWDIr    : (nTime x 1) measured EWDI trajectory. Values before t0 are
%              returned as NaN.

    n = size(sample, 1);
    EWDIr = NaN(n, 1);

    for i = t0+1:n
        Xt = sample(1:i, :);
        Xc = Xt - mean(Xt, 1);
        [~, ~, V] = svd(Xc, 'econ');
        loading_t = V(:, 1);
        EWDIr(i) = 1 - abs(loading_t' * loading0);
    end
end
