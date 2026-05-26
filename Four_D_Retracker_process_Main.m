%% ========================================================================
%  Run LBiST retracking and compare against OCOG / OCEAN / ICE / IT
%
%  Workflow:
%    1) Load precomputed 4D stack for one case study
%    2) Compute spatial and spatio-temporal masks
%    3) Reconstruct filtered waveforms
%    4) Apply retracking methods
%    5) Compute water levels and validation metrics
%    6) Run cycle-wise analysis
%    7) Save outputs
%% ========================================================================

clear; clc; close all;
warning off;

%% ----------------------- INITIALIZATION ----------------------------------
Lib_Initialization;

outlierThreshold = 0.8;
rangeResolution = 1280;
rangeResolutionRatio = rangeResolution / 128;

caseStudy = "LAKE MOHAVE AT DAVIS DAM, AZ-NV - 09422500";
load(sprintf('%d_%s.mat', rangeResolution, caseStudy));

%% ----------------------- BAYESIAN MASKS ----------------------------------
tic;

spatialTemporalMasks = Lib_Bayesian_PCA_FFT(Four_D_STACK); 



%% ----------------------- OUTPUT ARRAYS -----------------------------------
waterLevelSAMOSA = [];
waterLevelLBiST  = [];
waterLevelOCOG   = [];
waterLevelOCEAN  = [];
waterLevelICE    = [];
waterLevelInSitu = [];
waterLevelIT     = [];
timeFinal        = [];

wfSimulated = [];
wfOriginal  = [];

%% ----------------------- MAIN RETRACKING LOOP ----------------------------
for cycleIdx = 1:size(Four_D_STACK, 1)

    % Extract cycle data
    shiftBins       = Four_D_STACK{cycleIdx, 2};
    altL1BS         = Four_D_STACK{cycleIdx, 3};
    rangeL1BS       = Four_D_STACK{cycleIdx, 4};
    corrections     = Four_D_STACK{cycleIdx, 5};
    timeVec         = Four_D_STACK{cycleIdx, 6};
    cycleNumber     = Four_D_STACK{cycleIdx, 7}; 
    lats            = Four_D_STACK{cycleIdx, 8}; 
    lons            = Four_D_STACK{cycleIdx, 9}; 
    refWL           = Four_D_STACK{cycleIdx,10};
    ddmScaled       = Four_D_STACK{cycleIdx,11};
    rangeOCOG       = Four_D_STACK{cycleIdx,12};
    rangeOCEAN      = Four_D_STACK{cycleIdx,13};
    rangeICE        = Four_D_STACK{cycleIdx,18};

    nMeas = size(ddmScaled, 3);

    %% ------------------- BUILD L2-LIKE WAVEFORMS ------------------------
    wfCycle = zeros(nMeas, size(ddmScaled, 1));
    for kk = 1:nMeas
        wfCycle(kk, :) = nansum(ddmScaled(:, :, kk), 2)';
    end

    % Estimate common peak gate from waveform stack
    stackedWaveform = prod(wfCycle, 1, 'omitnan');
    [~, peakGate] = max(stackedWaveform);

    %% ------------------- PROCESS EACH MEASUREMENT -----------------------
    for kk = 1:nMeas

        peakGateThis = peakGate + shiftBins(kk);

        sliceOriginal = ddmScaled(:, :, kk);

        % Original waveform before masking
        wfOriginalNotShifted = nansum(sliceOriginal, 2);

        % Spatio-temporal mask collapsed over looks
        DDM_Simulated_NotShifted = sliceOriginal .* spatialTemporalMasks(:, :, cycleIdx);
        wfSimulated(kk, :) = nansum(DDM_Simulated_NotShifted, 2);

        % Normalize
        wfOriginalNotShifted = wfOriginalNotShifted ./ max(wfOriginalNotShifted );

        % Simulated / filtered waveform
        wfSimulated(kk, :) = wfSimulated(kk, :) ./ max(wfSimulated(kk, :) );
        wfSimulated(kk, :) = circshift(wfSimulated(kk, :), shiftBins(kk));

        % Original waveform with alignment
        wfOriginal(kk, :) = nansum(sliceOriginal, 2)';
        wfOriginal(kk, :) = circshift(wfOriginal(kk, :), shiftBins(kk));

        %% ---------------- RETRACK LBiST WAVEFORM -----------------------
        normalizedSimWF = wfSimulated(kk, :) ./ max(wfSimulated(kk, :) );

        thresh = 0.5;
        idx = find(normalizedSimWF < thresh & circshift(normalizedSimWF, -1) >= thresh, 1, 'first');

        if isempty(idx) || idx >= numel(normalizedSimWF)
            continue;
        end

        y1 = normalizedSimWF(idx);
        y2 = normalizedSimWF(idx + 1);
        frac = (thresh - y1) / (y2 - y1 );

        xHalfMax = idx + frac;
        xHalfMax = xHalfMax / rangeResolutionRatio;

        retrackCorrection = (xHalfMax - 44) * tau;
        rawElevation = altL1BS(kk) - rangeL1BS(kk);

        %% ---------------- OTHER RETRACKERS -----------------------------
        [rtrckrBinOCOG] = ocog_retracker(wfOriginal(kk, :) ./ max(wfOriginal(kk, :) ), 0.5, 0);
        [rtrckrBinIT]   = Retracker4_Thereshold(wfOriginal(kk, :) ./ max(wfOriginal(kk, :) ), 50, rangeResolution);

        rtrckrBinOCOG = rtrckrBinOCOG / rangeResolutionRatio;
        rtrckrBinIT   = rtrckrBinIT   / rangeResolutionRatio;

        if isnan(rtrckrBinIT)
            rtrckrBinIT = rtrckrBinOCOG;
        end

        % Downsampled waveform for SAMOSA+ if needed
        wfDown = downsample(wfOriginal(kk, :) ./ max(wfOriginal(kk, :) ), rangeResolutionRatio);

        % Placeholder: activate if SAMOSA+ retracker is available
        samosaPlus = 0;
        % samosaPlus = Lib_samosa_retracker(wfDown, rangeL1BS(kk), altL1BS(kk), ...
        %     corrections(kk), GEO, peakGateThis / rangeResolutionRatio);

        %% ---------------- WATER LEVELS --------------------------------
        waterLevelSAMOSA = [waterLevelSAMOSA samosaPlus];
        waterLevelLBiST  = [waterLevelLBiST rawElevation - retrackCorrection - corrections(kk)];
        waterLevelOCOG   = [waterLevelOCOG  rawElevation - ((rtrckrBinOCOG - 44) * tau) - corrections(kk)];
        waterLevelIT     = [waterLevelIT    rawElevation - ((rtrckrBinIT   - 44) * tau) - corrections(kk)];

        % Fill missing official retracker ranges
        if isnan(rangeOCOG(kk))
            rangeOCOG(kk) = rangeL1BS(kk) - ((rtrckrBinOCOG - 44) * tau);
        end

        if isnan(rangeICE(kk))
            rangeICE(kk) = rangeOCOG(kk);
        end

        if isnan(rangeOCEAN(kk))
            rangeOCEAN(kk) = rangeOCOG(kk);
        end

        waterLevelICE   = [waterLevelICE   altL1BS(kk) - rangeICE(kk)   - corrections(kk)];
        waterLevelOCEAN = [waterLevelOCEAN altL1BS(kk) - rangeOCEAN(kk) - corrections(kk)];

        waterLevelInSitu = [waterLevelInSitu refWL(kk)];
        timeFinal        = [timeFinal timeVec(kk)];
    end
end

toc;

%% ----------------------- BIAS CORRECTION ---------------------------------
biasLBiST  = median(waterLevelLBiST  - waterLevelInSitu);
biasOCOG   = median(waterLevelOCOG   - waterLevelInSitu);
biasSAMOSA = median(waterLevelSAMOSA - waterLevelInSitu);
biasOCEAN  = median(waterLevelOCEAN  - waterLevelInSitu);
biasIT     = median(waterLevelIT     - waterLevelInSitu);

% ICE uses valid subset only
nanIdx = isnan(waterLevelICE);
waterLevelICEValid = waterLevelICE(~nanIdx);
waterLevelInSituICE = waterLevelInSitu(~nanIdx);
biasICE = nanmedian(waterLevelICEValid - waterLevelInSituICE);

%% ----------------------- METRICS -----------------------------------------
metricSAMOSA = Metrics(waterLevelInSitu', waterLevelSAMOSA' - biasSAMOSA);
metricOCOG   = Metrics(waterLevelInSitu', waterLevelOCOG'   - biasOCOG);
metricLBiST  = Metrics(waterLevelInSitu', waterLevelLBiST'  - biasLBiST);
metricOCEAN  = Metrics(waterLevelOCEAN' - biasOCEAN, waterLevelInSitu');
metricICE    = Metrics(waterLevelICEValid' - biasICE, waterLevelInSituICE');
metricIT     = Metrics(waterLevelIT' - biasIT, waterLevelInSitu');

%% ----------------------- QUICK PLOT --------------------------------------
figure('Units', 'inches', 'Position', [1, 1, 10, 4]);
hold on;

plot(waterLevelLBiST - biasLBiST, 'b', 'LineWidth', 3);
plot(waterLevelOCOG  - biasOCOG,  'r', 'LineWidth', 2);
plot(waterLevelSAMOSA - biasSAMOSA, 'c', 'LineWidth', 2);
plot(waterLevelInSitu, 'k', 'LineWidth', 2);

title(caseStudy, 'Interpreter', 'none');
legend('LBiST', 'OCOG', 'SAMOSA', 'In-Situ', 'Location', 'best');
ylabel('WSH (m)');

tickIdx = 1:floor(numel(waterLevelLBiST)/10):numel(waterLevelLBiST);
xticks(tickIdx);
xticklabels(datestr(Tag_Time_AltBundle(timeFinal(tickIdx)), 'mmm yyyy'));
set(gca, 'FontSize', 16);
grid on;
set(gca, 'LooseInset', get(gca, 'TightInset'));

%% ----------------------- OUTLIER COUNTS ----------------------------------
nSamples = numel(timeFinal);

outliersLBiST = countOutliers(waterLevelLBiST, waterLevelInSitu, biasLBiST, outlierThreshold);
outliersOCOG  = countOutliers(waterLevelOCOG,  waterLevelInSitu, biasOCOG,  outlierThreshold);

%% ----------------------- PRINT RESULTS -----------------------------------
fprintf('************************\n');
fprintf('%s\n', caseStudy);
fprintf('range resolution: %d\n', rangeResolution);
fprintf('All meas 4D vs OCOG\n');
fprintf('        Correlation & %2.2f & %2.2f \\\\\n', metricLBiST.Corr, metricOCOG.Corr);
fprintf('        KGE & %2.2f & %2.2f \\\\\n', metricLBiST.KGE, metricOCOG.KGE);
fprintf('        NSE & %2.2f & %2.2f \\\\\n', metricLBiST.NSE, metricOCOG.NSE);
fprintf('        RMSE & %2.2f & %2.2f \\\\\n', metricLBiST.RMSE, metricOCOG.RMSE);
fprintf('        Number of outliers in %d samples & %d & %d \\\\\n', ...
    nSamples, outliersLBiST, outliersOCOG);

%% ----------------------- CYCLE ANALYSIS ----------------------------------
cycleAnalysis = true;
if cycleAnalysis
    Lib_Cycle_Analysis;
end

%% ----------------------- SAVE RESULTS ------------------------------------
save(sprintf('Final_Mean_Result_SP_%d_%s.mat', rangeResolution, caseStudy), ...
    'waterLevelLBiST', 'waterLevelICE', 'waterLevelOCEAN', ...
    'waterLevelOCOG', 'waterLevelSAMOSA', 'waterLevelInSitu', ...
    'waterLevelIT', 'biasLBiST', 'biasOCOG', 'biasSAMOSA', ...
    'biasOCEAN', 'biasICE', 'biasIT', ...
    'metricLBiST', 'metricOCOG', 'metricIT', 'metricSAMOSA', ...
    'metricOCEAN', 'metricICE', 'timeFinal', ...
    'Metric_OCOG_Cycle', 'Metric_LBiST_Cycle', 'Metric_SAMOSA_Cycle', ...
    'Metric_OCEAN_Cycle', 'Metric_ICE_Cycle', 'Metric_IT_Cycle');

%% ========================================================================
%  LOCAL FUNCTIONS
%% ========================================================================

function lut = readLUT(filename)
    fileID = fopen(filename, 'r');
    if fileID < 0
        error('Could not open LUT file: %s', filename);
    end
    lut = fscanf(fileID, '%f %f', [2 Inf])';
    fclose(fileID);
end

function nOutliers = countOutliers(waterLevel, waterLevelInSitu, bias, threshold)
    nOutliers = 0;
    for i = 1:numel(waterLevel)
        if abs(waterLevel(i) - (waterLevelInSitu(i) + bias)) > threshold
            nOutliers = nOutliers + 1;
        end
    end
end