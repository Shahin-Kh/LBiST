%% ========================================================================
%  Cycle-median analysis
%
%  Uses variable names consistent with the cleaned retracking script:
%    waterLevelLBiST, waterLevelOCOG, waterLevelSAMOSA, ...
%    waterLevelInSitu, timeFinal, caseStudy, outlierThreshold
%% ========================================================================

% Unique cycle timestamps
timesUnique = unique(timeFinal);
nCycles = numel(timesUnique);

% Preallocate cycle-level series
waterLevelCycle_LBiST  = nan(1, nCycles);
waterLevelCycle_OCOG   = nan(1, nCycles);
waterLevelCycle_SAMOSA = nan(1, nCycles);
waterLevelCycle_ICE    = nan(1, nCycles);
waterLevelCycle_IT     = nan(1, nCycles);
waterLevelCycle_OCEAN  = nan(1, nCycles);
waterLevelCycle_InSitu = nan(1, nCycles);
timeCycle              = nan(1, nCycles);

%% ----------------------- MEDIAN AGGREGATION ------------------------------
for i = 1:nCycles
    idx = (timeFinal == timesUnique(i));

    waterLevelCycle_LBiST(i)  = nanmedian(waterLevelLBiST(idx));
    waterLevelCycle_OCOG(i)   = nanmedian(waterLevelOCOG(idx));
    waterLevelCycle_InSitu(i) = nanmedian(waterLevelInSitu(idx));
    waterLevelCycle_SAMOSA(i) = nanmedian(waterLevelSAMOSA(idx));
    waterLevelCycle_ICE(i)    = nanmedian(waterLevelICE(idx));
    waterLevelCycle_OCEAN(i)  = nanmedian(waterLevelOCEAN(idx));
    waterLevelCycle_IT(i)     = nanmedian(waterLevelIT(idx));

    timeCycle(i) = timeFinal(find(idx, 1, 'first'));
end

%% ----------------------- BIAS ESTIMATION ---------------------------------
biasLBiST_Cycle  = nanmedian(waterLevelCycle_LBiST  - waterLevelCycle_InSitu);
biasOCOG_Cycle   = nanmedian(waterLevelCycle_OCOG   - waterLevelCycle_InSitu);
biasSAMOSA_Cycle = nanmedian(waterLevelCycle_SAMOSA - waterLevelCycle_InSitu);
biasICE_Cycle    = nanmedian(waterLevelCycle_ICE    - waterLevelCycle_InSitu);
biasOCEAN_Cycle  = nanmedian(waterLevelCycle_OCEAN  - waterLevelCycle_InSitu);
biasIT_Cycle     = nanmedian(waterLevelCycle_IT     - waterLevelCycle_InSitu);

%% ----------------------- OUTLIER COUNTS ----------------------------------
nCycleSamples = numel(waterLevelCycle_LBiST);

outliersLBiST_Cycle = countOutliers( ...
    waterLevelCycle_LBiST, waterLevelCycle_InSitu, biasLBiST_Cycle, outlierThreshold);

outliersOCOG_Cycle = countOutliers( ...
    waterLevelCycle_OCOG, waterLevelCycle_InSitu, biasOCOG_Cycle, outlierThreshold);

%% ----------------------- METRICS -----------------------------------------
Metric_OCOG_Cycle   = Metrics(waterLevelCycle_InSitu', waterLevelCycle_OCOG'   - biasOCOG_Cycle);
Metric_LBiST_Cycle  = Metrics(waterLevelCycle_InSitu', waterLevelCycle_LBiST'  - biasLBiST_Cycle);
Metric_SAMOSA_Cycle = Metrics(waterLevelCycle_InSitu', waterLevelCycle_SAMOSA' - biasSAMOSA_Cycle);
Metric_OCEAN_Cycle  = Metrics(waterLevelCycle_InSitu', waterLevelCycle_OCEAN'  - biasOCEAN_Cycle);
Metric_ICE_Cycle    = Metrics(waterLevelCycle_InSitu', waterLevelCycle_ICE'    - biasICE_Cycle);
Metric_IT_Cycle     = Metrics(waterLevelCycle_InSitu', waterLevelCycle_IT'     - biasIT_Cycle);

%% ----------------------- PLOT --------------------------------------------
figure('Units', 'inches', 'Position', [1, 1, 10, 4]);
hold on;

plot(waterLevelCycle_LBiST  - biasLBiST_Cycle,  'LineWidth', 2);
plot(waterLevelCycle_OCOG   - biasOCOG_Cycle,   'LineWidth', 2);
plot(waterLevelCycle_SAMOSA - biasSAMOSA_Cycle, 'LineWidth', 2);
plot(waterLevelCycle_InSitu, 'LineWidth', 2);

title(caseStudy, 'Interpreter', 'none');
legend('LBiST', 'OCOG', 'SAMOSA', 'In-Situ', 'Location', 'best');
ylabel('WSH (m)');

tickIdx = 1:floor(numel(waterLevelCycle_LBiST) / 10):numel(waterLevelCycle_LBiST);
xticks(tickIdx);

tickDates = Tag_Time_AltBundle(timeCycle(tickIdx));
xticklabels(datestr(tickDates, 'mmm yyyy'));

set(gca, 'FontSize', 16);
grid on;

%% ----------------------- PRINT RESULTS -----------------------------------
fprintf('%s\n', caseStudy);
fprintf('Cycle median 4D vs OCOG\n');
fprintf('        Correlation & %2.2f & %2.2f \\\\\n', Metric_LBiST_Cycle.Corr, Metric_OCOG_Cycle.Corr);
fprintf('        KGE & %2.2f & %2.2f \\\\\n', Metric_LBiST_Cycle.KGE, Metric_OCOG_Cycle.KGE);
fprintf('        NSE & %2.2f & %2.2f \\\\\n', Metric_LBiST_Cycle.NSE, Metric_OCOG_Cycle.NSE);
fprintf('        RMSE & %2.2f & %2.2f \\\\\n', Metric_LBiST_Cycle.RMSE, Metric_OCOG_Cycle.RMSE);
fprintf('        Number of outliers in %d samples & %d & %d \\\\\n', ...
    nCycleSamples, outliersLBiST_Cycle, outliersOCOG_Cycle);

clear timesUnique;

%% ========================================================================
%  Local function
%% ========================================================================
function nOutliers = countOutliers(waterLevel, waterLevelInSitu, bias, threshold)
    nOutliers = 0;
    for ii = 1:numel(waterLevel)
        if abs(waterLevel(ii) - (waterLevelInSitu(ii) + bias)) > threshold
            nOutliers = nOutliers + 1;
        end
    end
end