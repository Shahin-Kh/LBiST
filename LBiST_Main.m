%% LBiST Processing Script

clear;
close all;
warning off;

Lib_Initialization;

%% Configuration
cfg.generateDataSet     = false;
cfg.RunLBiST 	        = true;
cfg.binShiftOffset      = -10;

cfg.useSpatialTemporal  = true;
cfg.useSpatialOnly      = false;

cfg.useAllData          = true;
cfg.useCycleMedians     = false;
cfg.outlierThreshold    = 0.8;

cfg.highResolution      = true;   % 180 x 1280 pixels
cfg.debugFigures        = exist('debug_figures', 'var') && debug_figures;

%% Case study selection
% Select one or more case studies
caseStudies = "Missouri_River_06467500";

% Alternative examples:
% caseStudies = "Boca Res NR Truckee CA - 10344490";
% caseStudies = "LAKE MOHAVE AT DAVIS DAM, AZ-NV - 09422500";
% caseStudies = "Rathbun Lake near Rathbun, IA - 06903880";
% caseStudies = "Itaparica reservoir";
% caseStudies = "Lake Tahoe near Tahoe City CA - 10337000";
% caseStudies = "Lake Koshkonong Near Newville, WI - 05427235";
% caseStudies = "Columbia River at Stevenson, WA - 14128600";
% caseStudies = "Ohio River at Cannelton Dam at Cannelton, IN - 03303280";
% caseStudies = "Susquehanna River at Harrisburg, PA - 01570500";
% caseStudies = "Cheyenne River near Eagle Butte SD - 06439500";
% caseStudies = "COOSA River at Leesburg - 02399500";
% caseStudies = "Mohawk River at Vischer Ferry Dam NY 01356000";
% caseStudies = "Penobscot_River_01036390";
% caseStudies = "POTOMAC_River_01613000";
% caseStudies = "Upper_Red_Lake";
% caseStudies = "Lower_Red_Lake";

%% Main loop over case studies
for iCase = 1:numel(caseStudies)

    caseStudy = caseStudies(iCase);

    %% --------------------------------------------------------------------
    % Generate and save Four_D_STACK
    % ---------------------------------------------------------------------
    if cfg.generateDataSet

        Lib_Load_Datasets;

        x = 1;
        Four_D_STACK = cell(0, 18);

        for iCycle = 1:numel(cycles_L1BS)

            currentCycle = cycles_L1BS(iCycle);

            if ~ismember(currentCycle, cycles_alt_bundle) || ismember(currentCycle, outlier_cycles)
                continue;
            end

            ptrs = find(currentCycle == dataset.ObjVS.Gen.Sat.Lat.Hi.C.Cycle);

            % Keep unique latitude samples only
            lats_all = dataset.ObjVS.Gen.Sat.Lat.Hi.C.Signal(ptrs);
            [~, uniqueIdx] = unique(lats_all, 'stable');
            ptrs = ptrs(uniqueIdx);

            % Extract measurement-level variables
            lats       = dataset.ObjVS.Gen.Sat.Lat.Hi.C.Signal(ptrs);
            lons       = dataset.ObjVS.Gen.Sat.Lon.Hi.C.Signal(ptrs);
            xgm        = dataset.ObjVS.Gen.GeoH.XGM2019(ptrs);
            atm_dry    = dataset.ObjVS.Gen.Cor.Atm.DryTro.ECMWF.Hi.Signal(ptrs);
            atm_wet    = dataset.ObjVS.Gen.Cor.Atm.WetTro.ECMWF.Hi.Signal(ptrs);
            atm_ion    = dataset.ObjVS.Gen.Cor.Atm.Iono.Gim.Hi.Ku.Signal(ptrs);
            tar_ET     = dataset.ObjVS.Gen.Cor.Target.EarthTide.Hi.Signal(ptrs);
            tar_PT     = dataset.ObjVS.Gen.Cor.Target.PoleTide.Hi.Signal(ptrs);

            corr_tmp   = atm_dry + atm_wet + atm_ion + tar_ET + tar_PT + xgm;
            timeVec    = floor(datenum(Tag_Time_AltBundle(dataset.ObjVS.Raw.Mes.Wvf.power.SAR.Ku.Time(ptrs))));

            waveforms_L2          = dataset.ObjVS.Raw.Mes.Wvf.power.SAR.Ku.Signal(ptrs, :);
            altitude              = dataset.ObjVS.Raw.Sat.Alt.Hi.Ku.Signal(ptrs);
            rawElevation          = altitude - dataset.ObjVS.Raw.Mes.Rng.Tracker.SAR.Ku.Signal(ptrs);
            rawElevation_OCOG     = altitude - dataset.ObjVS.Raw.Mes.Rng.OCOG.SAR.Hi.Ku.Signal(ptrs);

            rangeTracker          = dataset.ObjVS.Raw.Mes.Rng.Tracker.SAR.Ku.Signal(ptrs);
            rangeOCOG             = dataset.ObjVS.Raw.Mes.Rng.OCOG.SAR.Hi.Ku.Signal(ptrs);
            rangeOCEAN            = dataset.ObjVS.Raw.Mes.Rng.Oce.SAR.Hi.Ku.Signal(ptrs);
            rangeICE              = dataset.ObjVS.Raw.Mes.Rng.Ice.SAR.Hi.Ku.Signal(ptrs);

            inputL1BS = string(L1BS(iCycle));
            lon_20Hz  = ncread(inputL1BS, 'lon_l1bs_echo_sar_ku');
            lat_20Hz  = ncread(inputL1BS, 'lat_l1bs_echo_sar_ku');

            % Bounding box around target region
            targetLatMin = min(lats) - 0.0132;
            targetLatMax = max(lats) + 0.0132;
            targetLonMin = min(lons) - 0.0132;
            targetLonMax = max(lons) + 0.0132;

            bboxMask = (lat_20Hz <= targetLatMax) & ...
                       (lat_20Hz >= targetLatMin) & ...
                       (lon_20Hz <= targetLonMax) & ...
                       (lon_20Hz >= targetLonMin);

            indBurst = find(bboxMask == 1);

            time_20Hz = ncread(inputL1BS, 'UTC_day_l1bs_echo_sar_ku', indBurst(1), 1);
            Time_Label = Tag_Time(time_20Hz); 

            if cfg.debugFigures
                figure;
                geoplot(lats, lons, 'r.-', 'MarkerSize', 30);
                geobasemap topographic;
            end

            ref_WL = nan(size(ptrs));
            DDM = [];
            DDM_with_scaling = [];
            wf = [];
            range_l1bs = nan(size(ptrs));
            alt_l1bs   = nan(size(ptrs));

            for kk = 1:numel(ptrs)

                [~, k] = min(CartesianDistance([lat_20Hz lon_20Hz], [lats(kk) lons(kk)]));

                ptrRef = find(timeVec(kk) == ref(:,4));
                if isempty(ptrRef)
                    continue;
                end
                ref_WL(kk) = mean(ref(ptrRef, 5));

                i_echo = ncread(inputL1BS, 'i_echoes_ku_l1bs_echo_sar_ku', [1 1 k], [128 180 1]);
                q_echo = ncread(inputL1BS, 'q_echoes_ku_l1bs_echo_sar_ku', [1 1 k], [128 180 1]);

                scale_factor = 1;
                DDM(:,:,kk) = ((scale_factor .* i_echo').^2 + (scale_factor .* q_echo').^2)';

                scale_factor = ncread(inputL1BS, 'iq_scale_factor_l1bs_echo_sar_ku', [1 k], [180 1]);
                DDM_with_scaling(:,:,kk) = ((scale_factor .* i_echo').^2 + (scale_factor .* q_echo').^2)';

                wf(kk,:) = nanmedian(DDM(:,:,kk), 2);

                if cfg.debugFigures
                    figure;
                    plot(wf(kk,:) / max(wf(kk,:)), 'r', 'LineWidth', 4);
                    pbaspect([3 1 1]);
                    axis off;
                    set(gca, 'XColor', 'none', 'YColor', 'none');
                    set(gcf, 'Color', 'none');
                    exportgraphics(gcf, sprintf('waveform_plot%d.png', kk), ...
                        'BackgroundColor', 'black', 'Resolution', 300);
                end

                range_l1bs(kk) = ncread(inputL1BS, 'range_ku_l1bs_echo_sar_ku', k, 1);
                alt_l1bs(kk)   = ncread(inputL1BS, 'alt_l1bs_echo_sar_ku', k, 1);
            end

            if isempty(ptrRef)
                continue;
            end

            % Store metadata
            Four_D_STACK{x,3}  = alt_l1bs;
            Four_D_STACK{x,4}  = range_l1bs;
            Four_D_STACK{x,5}  = corr_tmp;
            Four_D_STACK{x,6}  = timeVec;
            Four_D_STACK{x,7}  = currentCycle;
            Four_D_STACK{x,8}  = lats;
            Four_D_STACK{x,9}  = lons;
            Four_D_STACK{x,10} = ref_WL;

            rawElevation_cycle      = alt_l1bs - range_l1bs;
            rawElevation_cycle_mean = nanmedian(rawElevation_cycle);
            rawElevation_total_mean = nanmedian(Raw_elevation_total);

            bin_shift_global    = (rawElevation_cycle - rawElevation_total_mean) / tau;
            bin_shift_cyclewise = (rawElevation_cycle - rawElevation_cycle_mean) / tau; 

            % Upsample DDM from 128 to 1280 bins
            [~, ~, nMeas] = size(DDM);
            DDM_up = zeros(1280, 180, nMeas);
            DDM_with_scaling_up = zeros(1280, 180, nMeas);
            shift = zeros(1, numel(ptrs));

            for kk = 1:numel(ptrs)

                DDM_up(:,:,kk) = imresize(DDM(:,:,kk), [1280 180], 'bilinear');
                DDM_with_scaling_up(:,:,kk) = imresize(DDM_with_scaling(:,:,kk), [1280 180], 'bilinear');

                shiftValue = (bin_shift_global(kk) + cfg.binShiftOffset) * 10;

                if round(shiftValue) < 0
                    shift(kk) = -abs(round(shiftValue));
                else
                    shift(kk) = round(shiftValue);
                end

                DDM_up(:,:,kk) = circshift(DDM_up(:,:,kk), -shift(kk), 1);
                DDM_with_scaling_up(:,:,kk) = circshift(DDM_with_scaling_up(:,:,kk), -shift(kk), 1);
            end

            disp(iCycle);

            Four_D_STACK{x,1}  = DDM_up;
            Four_D_STACK{x,2}  = shift;
            Four_D_STACK{x,11} = DDM_with_scaling_up;
            Four_D_STACK{x,12} = rangeOCOG;
            Four_D_STACK{x,13} = rangeOCEAN;
            Four_D_STACK{x,14} = waveforms_L2;
            Four_D_STACK{x,15} = 0;
            Four_D_STACK{x,16} = 0;
            Four_D_STACK{x,17} = rawElevation_OCOG - corr_tmp;
            Four_D_STACK{x,18} = 0;

            clear alt_l1bs range_l1bs corr_tmp timeVec lats lons ref_WL ...
                  DDM shift DDM_with_scaling waveforms_L2 DDM_with_scaling_up DDM_up

            x = x + 1;
        end

        save(sprintf('%s.mat', caseStudy), 'Four_D_STACK', '-v7.3');
    end

    %% --------------------------------------------------------------------
    % Load and generate outputs / metrics / plots
    % ---------------------------------------------------------------------
    if cfg.RunLBiST

        load(sprintf('%s.mat', caseStudy));

        tic;
        [Spatial_masks, Spatial_Temporal_masks, ...
            reconstructed_Spatial_LB, reconstructed_Spatial_Temporal_LB] = ...
            Lib_Bayesian_PCA_FFT(Four_D_STACK); 

        Water_level_SAMOSA     = [];
        Water_level_new_method = [];
        Water_level_OCOG       = [];
        Water_level_OCEAN      = [];
        Water_level_ICE        = [];
        Water_level_InSitu     = [];
        Water_level_IT         = [];
        time_final             = [];
        pointers               = []; 
        MK = 1;

        for x = 1:size(Four_D_STACK, 1)

            DDM              = Four_D_STACK{x,1};
            shift            = Four_D_STACK{x,2} - (10 * cfg.binShiftOffset);
            Alt_l1bs         = Four_D_STACK{x,3};
            Range_l1bs       = Four_D_STACK{x,4};
            Corr_tmp         = Four_D_STACK{x,5};
            timeVec          = Four_D_STACK{x,6};
            cycleNumber      = Four_D_STACK{x,7}; 
            lats             = Four_D_STACK{x,8}; 
            lons             = Four_D_STACK{x,9}; 
            ref_WL           = Four_D_STACK{x,10};
            DDM_with_scaling = Four_D_STACK{x,11};

            Range_OCOG = Four_D_STACK{x,12};
            Range_OCEAN = Four_D_STACK{x,13};
            Range_ICE = Four_D_STACK{x,14};

            for kk = 1:size(DDM, 3)

                Stack_data = Four_D_STACK{x,1};
                Stack_data_original_with_scaling = Four_D_STACK{x,11};

                slice_original_with_scaling = Stack_data_original_with_scaling(:,:,kk);

                slice_original = Stack_data(:,:,kk);
                slice_original = slice_original / max(slice_original(:));

                masks = Spatial_Temporal_masks(:,:,x);
                masks = masks / max(masks(:));

                slice_modified = abs(masks .* slice_original);

                wf_simulated(kk,:) = nanmean(slice_modified, 2);
                wf_simulated(kk,:) = circshift(wf_simulated(kk,:), shift(kk));

                wf_original(kk,:) = nanmean(slice_original_with_scaling, 2)';
                wf_original(kk,:) = circshift(wf_original(kk,:), shift(kk));

                normalized_sim_wf = wf_simulated(kk,:) / max(wf_simulated(kk,:));
                thresh = 0.5;

                idx = find(normalized_sim_wf < thresh & ...
                           circshift(normalized_sim_wf, -1) >= thresh, ...
                           1, 'first');

                y1 = normalized_sim_wf(idx);
                y2 = normalized_sim_wf(idx + 1);

                frac = (thresh - y1) / (y2 - y1);
                x_half_max = idx + frac;
                x_half_max = x_half_max / 10;

                retrack_correction = (x_half_max - 44) * tau;
                Raw_Elevation = Alt_l1bs(kk) - Range_l1bs(kk);

                [rtrckr_bin_OCOG] = ocog_retracker(wf_original(kk,:) / max(wf_original(kk,:)), 0.5, 0);
                [rtrckr_bin_IT]   = Retracker4_Thereshold(wf_original(kk,:) / max(wf_original(kk,:)), 50, 128);

                rtrckr_bin_OCOG = rtrckr_bin_OCOG / 10;
                rtrckr_bin_IT   = rtrckr_bin_IT / 10;

                % SAMOSA retracker
                GEO.LAT = 0;
                GEO.LON = 0;
                GEO.Height = Alt_l1bs(kk);
                GEO.Vs = 0;
                GEO.Hrate = 0;
                GEO.Pitch = 0;
                GEO.Roll = 0;
                GEO.nu = 0;
                GEO.track_sign = 0;

                wf_down = downsample(wf_original(kk,:), 10);
                samosa_plus = Lib_samosa_retracker( ...
                    wf_down / max(wf_down), ...
                    Range_l1bs(kk), ...
                    Alt_l1bs(kk), ...
                    Corr_tmp(kk), ...
                    GEO);

                Water_level_SAMOSA = [Water_level_SAMOSA samosa_plus]; 

                Water_level_new_method = [Water_level_new_method ...
                    Raw_Elevation - retrack_correction - Corr_tmp(kk)]; 

                Water_level_OCOG = [Water_level_OCOG ...
                    Raw_Elevation - ((rtrckr_bin_OCOG - 44) * tau) - Corr_tmp(kk)]; 

                % Preserved exactly as in original script
                Water_level_IT = [Water_level_OCOG ...
                    Raw_Elevation - ((rtrckr_bin_IT - 44) * tau) - Corr_tmp(kk)]; 

                Water_level_OCEAN = [Water_level_OCEAN ...
                    Alt_l1bs(kk) - Range_OCEAN(kk) - Corr_tmp(kk)]; 

                Water_level_ICE = [Water_level_ICE ...
                    Alt_l1bs(kk) - Range_ICE(kk) - Corr_tmp(kk)]; 

                Water_level_InSitu = [Water_level_InSitu ref_WL(kk)]; 
                time_final = [time_final timeVec(kk)]; 

                MK = MK + 1;
            end
        end
        toc;

        %% Metrics and summary
        if cfg.useAllData

            bias  = median(Water_level_new_method - Water_level_InSitu);
            bias2 = median(Water_level_OCOG - Water_level_InSitu);
            bias3 = median(Water_level_SAMOSA - Water_level_InSitu);

            Metric_SAMOSA = Metrics(Water_level_InSitu', Water_level_SAMOSA' - bias3);
            Metric_OCOG   = Metrics(Water_level_InSitu', Water_level_OCOG' - bias2);
            Metric_LBiST = Metrics(Water_level_InSitu', Water_level_new_method' - bias);

            Metric_OCEAN  = Metrics(Water_level_OCEAN' - bias2, Water_level_InSitu');
            Metric_ICE    = Metrics(Water_level_ICE' - bias2, Water_level_InSitu');

            save(strcat('Lake_', caseStudy), ...
                'Water_level_new_method', 'Water_level_OCOG', 'Water_level_SAMOSA', ...
                'Water_level_InSitu', 'bias', 'bias2', 'bias3', ...
                'Metric_LBiST', 'Metric_OCOG', 'Metric_SAMOSA', 'time_final');

            %% Plot time series
            figure('Units', 'inches', 'Position', [1, 1, 10, 4]);
            hold on;

            plot(Water_level_new_method - bias, 'Color', 'b', 'LineWidth', 3);
            plot(Water_level_OCOG - bias2, 'Color', 'r', 'LineWidth', 2);
            plot(Water_level_SAMOSA - bias3, 'Color', 'c', 'LineWidth', 2);
            plot(Water_level_InSitu, 'Color', 'k', 'LineWidth', 2);

            title(caseStudy);
            legend('4D', 'OCOG', 'SAMOSA', 'In-Situ');
            ylabel('WSH (m)');

            tmpTicks = 1:floor(numel(Water_level_new_method) / 10):numel(Water_level_new_method);
            xticks(tmpTicks);

            xlb = Tag_Time_AltBundle(time_final(tmpTicks));
            set(gca, 'FontSize', 16);
            xticklabels(datestr(xlb, 'mmm yyyy'));
            grid on;

            set(gca, 'LooseInset', get(gca, 'TightInset'));

            %% Outlier counts
            nSamples = numel(time_final);

            outliersNew = find(abs(Water_level_new_method - (Water_level_InSitu + bias)) > cfg.outlierThreshold);
            outliersOCOG = find(abs(Water_level_OCOG - (Water_level_InSitu + bias2)) > cfg.outlierThreshold);

            Outliers_New_Method = numel(outliersNew);
            Outliers_OCOG = numel(outliersOCOG);

            Metric_OCOG   = Metrics(Water_level_InSitu', Water_level_OCOG' - bias2);
            Metric_LBiST = Metrics(Water_level_InSitu', Water_level_new_method' - bias);

            fprintf('%s\n', caseStudy);
            fprintf('All meas 4D vs OCOG\n');
            fprintf('        Correlation & %2.2f & %2.2f \\\\ \n', Metric_LBiST.Corr, Metric_OCOG.Corr);
            fprintf('        KGE & %2.2f & %2.2f \\\\ \n', Metric_LBiST.KGE, Metric_OCOG.KGE);
            fprintf('        NSE & %2.2f & %2.2f \\\\ \n', Metric_LBiST.NSE, Metric_OCOG.NSE);
            fprintf('        RMSE & %2.2f & %2.2f \\\\ \n', Metric_LBiST.RMSE, Metric_OCOG.RMSE);
            fprintf('        Number of outliers in %d samples & %d & %d  \\\\ \n', ...
                nSamples, Outliers_New_Method, Outliers_OCOG);

            %% Cycle-wise analysis
            Cycle_Analysis = true;
            if Cycle_Analysis
                Lib_Cycle_Analysis2;
            end

            %% Optional print blocks
            printMetrics = false;
            if printMetrics
                fprintf('%s', caseStudy);
                fprintf('************************\n');
                fprintf('4D\n');
                fprintf('Correlation: %2.2f\n', Metric_LBiST.Corr);
                fprintf('KGE: %2.2f\n', Metric_LBiST.KGE);
                fprintf('NSE: %2.2f\n', Metric_LBiST.NSE);
                fprintf('RMSE: %2.2f\n', Metric_LBiST.RMSE);
                fprintf('number of outliers: %d in %d samples\n', Outliers_New_Method, nSamples);

                fprintf('************************\n');
                fprintf('OCOG Retracker\n');
                fprintf('Correlation: %2.2f\n', Metric_OCOG.Corr);
                fprintf('KGE: %2.2f\n', Metric_OCOG.KGE);
                fprintf('NSE: %2.2f\n', Metric_OCOG.NSE);
                fprintf('RMSE: %2.2f\n', Metric_OCOG.RMSE);
                fprintf('number of outliers: %d in %d samples\n', Outliers_OCOG, nSamples);
            end

            exportLatex = false;
            if exportLatex
                fprintf('%s', caseStudy);
                fprintf('        Correlation & %2.2f & %2.2f \\\\ \n', Metric_LBiST.Corr, Metric_OCOG.Corr);
                fprintf('        KGE & %2.2f & %2.2f \\\\ \n', Metric_LBiST.KGE, Metric_OCOG.KGE);
                fprintf('        NSE & %2.2f & %2.2f \\\\ \n', Metric_LBiST.NSE, Metric_OCOG.NSE);
                fprintf('        RMSE & %2.2f & %2.2f \\\\ \n', Metric_LBiST.RMSE, Metric_OCOG.RMSE);
                fprintf('        Number of outliers in %d samples & %d & %d  \\\\ \n', ...
                    nSamples, Outliers_New_Method, Outliers_OCOG);
            end
        end
    end
end