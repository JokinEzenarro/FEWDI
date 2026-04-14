% =========================================================================
% ewdi_reference.m
%
% Author: Jokin Ezenarro
% Contact: jokin@food.ku.dk
% Last modified: 14/4/2026
%
% License/usage:
% Free to use and modify for non-commercial purposes, provided that proper
% citation is given to the original work and repository.
% =========================================================================

function loading0 = ewdi_reference(Xref)
% EWDI_REFERENCE Build the reference loading (PC1) from reference spectra.
%
% INPUT
%   Xref : (nSamples x nVariables) matrix containing the reference spectra.
%
% OUTPUT
%   loading0 : first loading vector obtained from SVD-based PCA on mean-
%              centered reference spectra.

    Xc = Xref - mean(Xref, 1);
    [~, ~, V] = svd(Xc, 'econ');
    loading0 = V(:, 1);
end
