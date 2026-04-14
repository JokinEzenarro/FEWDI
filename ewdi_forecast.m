% =========================================================================
% ewdi_forecast.m
%
% Author: Jokin Ezenarro
% Contact: jokin@food.ku.dk
% Last modified: 14/4/2026
%
% License/usage:
% Free to use and modify for non-commercial purposes, provided that proper
% citation is given to the original work and repository.
% =========================================================================

function results = ewdi_forecast(sample, loading0, times, params)
% EWDI_FORECAST Forecast EWDI with PCA score extrapolation and Monte Carlo
% uncertainty propagation.
%
% INPUTS
%   sample   : (nTime x nVariables) matrix with one spectral trajectory
%   loading0 : reference loading vector
%   times    : (nTime x 1) numeric time vector in hours
%   params   : structure with fields
%       .t0           : minimum number of points before forecasting
%       .trainWindow  : calibration window width in hours
%       .moveStep     : step in hours between forecasting origins
%       .predN        : number of forecasted points at each origin
%       .predStep     : spacing in hours between future predictions
%       .nPC          : number of PCs used for score forecasting
%       .polyOrder    : polynomial order for score extrapolation
%       .kCI          : uncertainty band multiplier, e.g. 2
%       .nMC          : number of Monte Carlo draws
%
% OUTPUT
%   results : structure containing
%       .EWDIr        measured EWDI trajectory
%       .predMean     predicted EWDI means by origin/horizon
%       .predSD       predicted EWDI standard deviations by origin/horizon
%       .predTimes    future time grids by origin
%       .originTimes  forecast origin times
%       .pc1Mean      predicted PC1 mean
%       .pc1SD        predicted PC1 standard deviation

    if size(sample, 1) ~= numel(times)
        error('Number of rows in sample must match length of times.');
    end
    times = times(:);

    EWDIr = ewdi_real(sample, loading0, params.t0);

    tStartForecast = times(params.t0 + 1);
    originTimes = (tStartForecast : params.moveStep : times(end)).';

    nOrigins = numel(originTimes);
    predMean = NaN(params.predN, nOrigins);
    predSD   = NaN(params.predN, nOrigins);
    predTimes = NaN(params.predN, nOrigins);
    pc1Mean = NaN(params.predN, nOrigins);
    pc1SD   = NaN(params.predN, nOrigins);

    rng(1);

    for k = 1:nOrigins
        tNow = originTimes(k);
        idxW = find(times >= (tNow - params.trainWindow) & times <= tNow);

        minPts = max(params.nPC + 1, params.polyOrder + 3);
        if numel(idxW) < minPts
            continue;
        end

        Xt = sample(idxW, :);
        tw = times(idxW);
        iEnd = idxW(end);

        tNew = tNow + params.predStep * (1:params.predN).';
        predTimes(:, k) = tNew;

        mW = mean(Xt, 1);
        XcW = Xt - mW;
        [U, S, V] = svd(XcW, 'econ');
        T = U * S;
        P = V(:, 1:params.nPC);

        Xreg = local_poly_design(tw, params.polyOrder);
        Y = T(:, 1:params.nPC);
        B = Xreg \ Y;

        Xnew = local_poly_design(tNew, params.polyOrder);
        Tpred = Xnew * B;

        [Sigma_eps, hStar] = local_score_uncertainty(Xreg, Y, Xnew);

        muPC1 = Tpred(:, 1);
        sdPC1 = zeros(params.predN, 1);
        for j = 1:params.predN
            Sigma_pred_j = (1 + hStar(j)) * Sigma_eps;
            Sigma_pred_j = (Sigma_pred_j + Sigma_pred_j') / 2;
            sdPC1(j) = sqrt(max(Sigma_pred_j(1,1), 0));
        end
        pc1Mean(:, k) = muPC1;
        pc1SD(:, k) = sdPC1;

        Lcell = cell(params.predN, 1);
        for j = 1:params.predN
            Sigma_pred_j = (1 + hStar(j)) * Sigma_eps;
            Sigma_pred_j = (Sigma_pred_j + Sigma_pred_j') / 2;
            [L, pchol] = chol(Sigma_pred_j, 'lower');
            if pchol ~= 0
                jitter = 1e-10 * max(trace(Sigma_pred_j) / params.nPC, 1);
                L = chol(Sigma_pred_j + jitter * eye(params.nPC), 'lower');
            end
            Lcell{j} = L;
        end

        for j = 1:params.predN
            simEWDI = NaN(params.nMC, 1);

            for s = 1:params.nMC
                Tsim_path = zeros(j, params.nPC);
                for jj = 1:j
                    z = randn(1, params.nPC);
                    Tsim_path(jj, :) = Tpred(jj, 1:params.nPC) + z * Lcell{jj}.';
                end

                Xsim_path = mW + Tsim_path * P.';
                Xmc = [sample(1:iEnd, :); Xsim_path];
                XcMC = Xmc - mean(Xmc, 1);
                [~, ~, Vmc] = svd(XcMC, 'econ');
                loading_t = Vmc(:, 1);
                simEWDI(s) = 1 - abs(loading_t' * loading0);
            end

            predMean(j, k) = mean(simEWDI, 'omitnan');
            predSD(j, k) = std(simEWDI, 0, 'omitnan');
        end
    end

    results = struct();
    results.EWDIr = EWDIr;
    results.predMean = predMean;
    results.predSD = predSD;
    results.predTimes = predTimes;
    results.originTimes = originTimes;
    results.pc1Mean = pc1Mean;
    results.pc1SD = pc1SD;
end

function X = local_poly_design(t, polyOrder)
    t = t(:);
    X = zeros(numel(t), polyOrder + 1);
    for p = 0:polyOrder
        X(:, p+1) = t.^p;
    end
end

function [Sigma_eps, hStar] = local_score_uncertainty(Xreg, Y, Xnew)
    B = Xreg \ Y;
    Yhat = Xreg * B;
    R = Y - Yhat;
    dof = max(size(Y,1) - size(Xreg,2), 1);
    Sigma_eps = (R' * R) / dof;
    XtX_inv = inv(Xreg' * Xreg);
    hStar = sum((Xnew * XtX_inv) .* Xnew, 2);
end
