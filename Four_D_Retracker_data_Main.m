%% ========================================================================
%  Build 4D LBiST stack from Sentinel-3 L1BS data
%
%  This script:
%    1) Loads a case study and associated datasets
%    2) Extracts 20 Hz measurements for valid cycles
%    3) Matches L2 measurements to L1BS echoes
%    4) Builds look-bin maps (DDMs)
%    5) Upsamples and aligns them in range
%    6) Stores all required variables in Four_D_STACK
%
%% ========================================================================

clear; clc; close all;
warning off;

%% ----------------------- USER SETTINGS -----------------------------------
rangeResolution = 1280;
nativeRangeBins = 128;
rangeResolutionRatio = rangeResolution / nativeRangeBins;

CaseStudy = "LAKE MOHAVE AT DAVIS DAM, AZ-NV - 09422500";

%% ----------------------- INITIALIZATION ----------------------------------
Lib_Initialization;
Lib_Load_Datasets;

% Optional:
% plot_case_study(CaseStudy, dataset, USGS_POS);

%% ----------------------- OUTPUT CONTAINER --------------------------------
Four_D_STACK = {};
stackIndex = 1;

%% ----------------------- MAIN LOOP ---------------------------------------
for cycleIdx = 1:numel(cycles_L1BS)

    fprintf('Processing cycle %d / %d\n', cycleIdx, numel(cycles_L1BS));

    currentCycle = cycles_L1BS(cycleIdx);

    % Skip unavailable / invalid cycles
    if ~ismember(currentCycle, cycles_alt_bundle)
        continue;
    end

    % Find all 20 Hz measurements belonging to current cycle
    ptrs = find(currentCycle == dataset.ObjVS.Raw.Sat.Lat.Hi.Ku.Cycle);
    if isempty(ptrs)
        continue;
    end

    %% -------------------- EXTRACT L2 / AUX DATA --------------------------
    lats = dataset.ObjVS.Raw.Sat.Lat.Hi.Ku.Signal(ptrs);
    lons = dataset.ObjVS.Raw.Sat.Lon.Hi.Ku.Signal(ptrs);

    xgm      = dataset.ObjVS.Gen.GeoH.XGM2019(ptrs);
    atmDry   = dataset.ObjVS.Gen.Cor.Atm.DryTro.ECMWF.Hi.Signal(ptrs);
    atmWet   = dataset.ObjVS.Gen.Cor.Atm.WetTro.ECMWF.Hi.Signal(ptrs);
    atmIono  = dataset.ObjVS.Gen.Cor.Atm.Iono.Gim.Hi.Ku.Signal(ptrs);

    earthTide = dataset.ObjVS.Gen.Cor.Target.EarthTide.Hi.Signal(ptrs);
    poleTide  = dataset.ObjVS.Gen.Cor.Target.PoleTide.Hi.Signal(ptrs);

    corrections = atmDry + atmWet + atmIono + earthTide + poleTide + xgm;

    timeVec = floor(datenum(Tag_Time_AltBundle( ...
        dataset.ObjVS.Raw.Mes.Wvf.power.SAR.Ku.Time(ptrs))));

    waveformsL2 = dataset.ObjVS.Raw.Mes.Wvf.power.SAR.Ku.Signal(ptrs, :);
    altitude    = dataset.ObjVS.Raw.Sat.Alt.Hi.Ku.Signal(ptrs);

    trackerRange = dataset.ObjVS.Raw.Mes.Rng.Tracker.SAR.Ku.Signal(ptrs);
    rangeOCOG    = dataset.ObjVS.Raw.Mes.Rng.OCOG.SAR.Hi.Ku.Signal(ptrs);
    rangeOCEAN   = dataset.ObjVS.Raw.Mes.Rng.Oce.SAR.Hi.Ku.Signal(ptrs);
    rangeICE     = dataset.ObjVS.Raw.Mes.Rng.Ice.SAR.Hi.Ku.Signal(ptrs);

    rawElevation      = altitude - trackerRange;
    rawElevationOCOG  = altitude - rangeOCOG;

    %% -------------------- READ L1BS FILE ---------------------------------
    inputL1BS = string(L1BS(cycleIdx));

    lon20Hz = ncread(inputL1BS, 'lon_l1bs_echo_sar_ku');
    lat20Hz = ncread(inputL1BS, 'lat_l1bs_echo_sar_ku');

    % Bounding box around target measurements
    bboxMargin = 0.0132;
    targetLatMin = min(lats) - bboxMargin;
    targetLatMax = max(lats) + bboxMargin;
    targetLonMin = min(lons) - bboxMargin;
    targetLonMax = max(lons) + bboxMargin;

    inBBox = (lat20Hz <= targetLatMax) & ...
             (lat20Hz >= targetLatMin) & ...
             (lon20Hz <= targetLonMax) & ...
             (lon20Hz >= targetLonMin);

    burstIdx = find(inBBox);
    if isempty(burstIdx)
        continue;
    end

    % Optional metadata
    time20Hz = ncread(inputL1BS, 'UTC_day_l1bs_echo_sar_ku', burstIdx(1), 1);
    Time_Label = Tag_Time(time20Hz);

    %% -------------------- PREALLOCATE CYCLE ARRAYS -----------------------
    nMeas = numel(ptrs);

    refWL         = nan(nMeas, 1);
    rangeL1BS     = nan(nMeas, 1);
    altL1BS       = nan(nMeas, 1);
    ddmRaw        = nan(nativeRangeBins, 180, nMeas);
    ddmScaled     = nan(nativeRangeBins, 180, nMeas);

    %% -------------------- BUILD DDMs -------------------------------------
    for kk = 1:nMeas

        % Find nearest L1BS echo to current L2 point
        [~, nearestIdx] = min(CartesianDistance([lat20Hz lon20Hz], [lats(kk) lons(kk)]));

        agc = ncread(inputL1BS, 'agc_ku_l1bs_echo_plrm', nearestIdx, 1);

        % Reference WL lookup
        refPtr = find(timeVec(kk) == ref(:, 4));
        if isempty(refPtr)
            continue;
        end
        refWL(kk) = mean(ref(refPtr, 5));

        % Read I/Q echoes
        iEcho = ncread(inputL1BS, 'i_echoes_ku_l1bs_echo_sar_ku', [1 1 nearestIdx], [128 180 1]);
        qEcho = ncread(inputL1BS, 'q_echoes_ku_l1bs_echo_sar_ku', [1 1 nearestIdx], [128 180 1]);

        % Unscaled power
        ddmRaw(:, :, kk) = (iEcho'.^2 + qEcho'.^2)';

        % Scaled power
        iqScale = ncread(inputL1BS, 'iq_scale_factor_l1bs_echo_sar_ku', [1 nearestIdx], [180 1]);
        ddmScaled(:, :, kk) = ((iqScale .* iEcho').^2 + (iqScale .* qEcho').^2)';
        ddmScaled(:, :, kk) = ddmScaled(:, :, kk) ./ 10^(agc / 10);

        % Geophysical range / altitude
        rangeL1BS(kk) = ncread(inputL1BS, 'range_ku_l1bs_echo_sar_ku', nearestIdx, 1);
        altL1BS(kk)   = ncread(inputL1BS, 'alt_l1bs_echo_sar_ku', nearestIdx, 1);
    end

    % Skip cycle if all reference values are missing
    if all(isnan(refWL))
        continue;
    end

    %% -------------------- RANGE ALIGNMENT --------------------------------
    rawElevationCycle      = altL1BS - rangeL1BS;
    rawElevationGlobalMean = nanmedian(Raw_elevation_total);  

    binShiftGlobal = (rawElevationCycle - rawElevationGlobalMean) / tau;

    ddmRawUp    = nan(rangeResolution, 180, nMeas);
    ddmScaledUp = nan(rangeResolution, 180, nMeas);
    shiftBins   = zeros(nMeas, 1);

    for kk = 1:nMeas
        % Upsample in range direction
        ddmRawUp(:, :, kk)    = imresize(ddmRaw(:, :, kk), [rangeResolution 180], 'bilinear');
        ddmScaledUp(:, :, kk) = imresize(ddmScaled(:, :, kk), [rangeResolution 180], 'bilinear');

        % Convert shift to upsampled grid
        shiftVal = binShiftGlobal(kk) * rangeResolutionRatio;
        shiftBins(kk) = round(shiftVal);

        % Align
        ddmRawUp(:, :, kk)    = circshift(ddmRawUp(:, :, kk),    -shiftBins(kk), 1);
        ddmScaledUp(:, :, kk) = circshift(ddmScaledUp(:, :, kk), -shiftBins(kk), 1);
    end

    %% -------------------- STORE RESULTS ----------------------------------
    Four_D_STACK{stackIndex, 1}  = 0;                  % placeholder for raw DDM upsampled
    Four_D_STACK{stackIndex, 2}  = shiftBins;
    Four_D_STACK{stackIndex, 3}  = altL1BS;
    Four_D_STACK{stackIndex, 4}  = rangeL1BS;
    Four_D_STACK{stackIndex, 5}  = corrections;
    Four_D_STACK{stackIndex, 6}  = timeVec;
    Four_D_STACK{stackIndex, 7}  = currentCycle;
    Four_D_STACK{stackIndex, 8}  = lats;
    Four_D_STACK{stackIndex, 9}  = lons;
    Four_D_STACK{stackIndex,10}  = refWL;
    Four_D_STACK{stackIndex,11}  = ddmScaledUp;
    Four_D_STACK{stackIndex,12}  = rangeOCOG;
    Four_D_STACK{stackIndex,13}  = rangeOCEAN;
    Four_D_STACK{stackIndex,14}  = waveformsL2;
    Four_D_STACK{stackIndex,15}  = 0;                  % placeholder
    Four_D_STACK{stackIndex,16}  = 0;                  % placeholder
    Four_D_STACK{stackIndex,17}  = rawElevationOCOG - corrections;
    Four_D_STACK{stackIndex,18}  = rangeICE;

    stackIndex = stackIndex + 1;

end

%% ----------------------- SAVE ---------------------------------
save(sprintf('%d_%s.mat', rangeResolution, CaseStudy), 'Four_D_STACK', '-v7.3');

