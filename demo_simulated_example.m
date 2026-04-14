% =========================================================================
% demo_paper_example.m
%
% Author: Jokin Ezenarro
% Contact: jokin@food.ku.dk
% Last modified: 14/4/2026
%
% License/usage:
% Free to use and modify for non-commercial purposes, provided that proper
% citation is given to the original work and repository.
% =========================================================================

% DEMO_PAPER_EXAMPLE
% Live demonstration of EWDI forecast evolution with shaded uncertainty.

clear; clc; close all;

% -------------------------------------------------------------------------
% Synthetic MIR-like fermentation data
% -------------------------------------------------------------------------
times = (0:4:208).';
nT = numel(times);
wn = linspace(4000, 850, 600);
sample = zeros(nT, numel(wn));

rng(2);
for i = 1:nT
    t = times(i);
    prog = 1 / (1 + exp(-(t - 80)/18));

    baseline = 0.15 + 0.02 * sin(wn/420) + 0.01 * cos(wn/170);
    band1 = 0.55 * exp(-0.5 * ((wn - 3300)/170).^2);
    band2 = (0.08 + 0.22*prog) * exp(-0.5 * ((wn - 2970)/55).^2);
    band3 = (0.05 + 0.16*prog) * exp(-0.5 * ((wn - 2885)/45).^2);
    band4 = (0.40 - 0.24*prog) * exp(-0.5 * ((wn - 1085)/42).^2);
    band5 = (0.22 - 0.11*prog) * exp(-0.5 * ((wn - 1045)/28).^2);
    band6 = (0.11 + 0.08*prog) * exp(-0.5 * ((wn - 1720)/35).^2);
    band7 = (0.10 + 0.03*sin(t/35)) * exp(-0.5 * ((wn - 1450)/30).^2);
    band8 = (0.06 + 0.05*prog) * exp(-0.5 * ((wn - 1260)/38).^2);
    drift = 0.015 * (t / max(times)) * sin(wn/95);
    noise = 0.003 * randn(1, numel(wn));

    sample(i, :) = baseline + band1 + band2 + band3 + band4 + band5 + ...
                   band6 + band7 + band8 + drift + noise;
end

% Reference and measured EWDI
Xref = sample(1:12, :);
loading0 = ewdi_reference(Xref);

t0 = 8;
params = struct();
params.t0 = t0;
params.trainWindow = 36;
params.moveStep = 4;
params.predN = 6;
params.predStep = 4;
params.nPC = 3;
params.polyOrder = 2;
params.kCI = 2;
params.nMC = 200;

results = ewdi_forecast(sample, loading0, times, params);

% -------------------------------------------------------------------------
% Static figure of synthetic MIR-like data
% -------------------------------------------------------------------------
figure('Color', 'w', 'Position', [50 80 900 350]);
imagesc(wn, times, sample);
set(gca, 'XDir', 'reverse');
xlabel('Wavenumber (cm^{-1})');
ylabel('Time (h)');
title('Synthetic MIR-like spectral evolution during wine fermentation');
colorbar;

% -------------------------------------------------------------------------
% Live EWDI forecast evolution with shaded uncertainty
% -------------------------------------------------------------------------
figure('Color', 'w', 'Position', [100 100 900 500]);
hold on; box on;
plot(times, results.EWDIr, 'k.-', 'LineWidth', 1.1, 'DisplayName', 'Measured EWDI');
xlabel('Time (h)');
ylabel('EWDI');
title('Live EWDI forecast evolution with uncertainty');
grid on;
ylim([0 max(0.8, nanmax(results.EWDIr) + 0.2)]);

hWindow = plot(NaN, NaN, 'bo-', 'LineWidth', 1.0, 'DisplayName', 'Calibration window');
hPred = plot(NaN, NaN, 'r.-', 'LineWidth', 1.4, 'DisplayName', 'Forecast mean');
hBand = fill(NaN, NaN, [1 0.6 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.35, ...
             'DisplayName', 'Forecast uncertainty');
hNow = xline(NaN, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.0, ...
             'DisplayName', 'Forecast origin');
legend('Location', 'northwest');

originTimes = results.originTimes;
for k = 1:numel(originTimes)
    valid = ~isnan(results.predMean(:, k)) & ~isnan(results.predTimes(:, k));
    if ~any(valid)
        continue;
    end

    tNow = originTimes(k);
    idxW = find(times >= (tNow - params.trainWindow) & times <= tNow);

    mu = results.predMean(valid, k);
    sd = results.predSD(valid, k);
    tt = results.predTimes(valid, k);

    set(hWindow, 'XData', times(idxW), 'YData', results.EWDIr(idxW));
    set(hPred, 'XData', tt, 'YData', mu);
    set(hBand, 'XData', [tt; flipud(tt)], ...
               'YData', [mu - params.kCI*sd; flipud(mu + params.kCI*sd)]);
    hNow.Value = tNow;

    drawnow;
    pause(0.20);
end
