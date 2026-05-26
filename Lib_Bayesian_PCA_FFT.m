function spatialTemporalMasks= ...
          Lib_Bayesian_PCA_FFT(ddmStack)
%% ========================================================================
%  Lib_Bayesian_PCA_FFT
%
%  Builds spatial and spatio-temporal masks from a 4D LBiST stack.
%
%  Inputs
%  ------
%  ddmStack : cell array
%      Cell array containing cycle-wise look-bin maps in ddmStack{cycle,11}
%
%  Outputs
%  -------
%  spatialTemporalMasks          : [range x look x cycle]
%
%  Method
%  ------
%  1) For each cycle, extract the dominant spatial structure via PCA
%  2) Stack reconstructed spatial maps across cycles
%  3) Estimate seasonal/annual temporal behavior with FFT
%  4) Combine seasonal and spatial patterns into a spatio-temporal mask
%% ========================================================================

    %% -------------------- INITIAL DIMENSIONS ----------------------------
    numCycles = size(ddmStack, 1);
    [numRows, numCols, ~] = size(ddmStack{1, 11});

    reconstructedSpatialLB         = zeros(numRows, numCols, numCycles);
    spatialMasks                   = zeros(numRows, numCols, numCycles);
    spatialTemporalMasks           = zeros(numRows, numCols, numCycles);
    reconstructedSpatialTemporalLB = zeros(numRows, numCols, numCycles);

    %% -------------------- SPATIAL PCA PER CYCLE -------------------------
    for cycle = 1:numCycles

        currentDDM = double(ddmStack{cycle, 11});
        [numRows, numCols, numMeas] = size(currentDDM);

        filteredDDM = currentDDM;

        % Reshape to [pixels x measurements]
        dataMatrix = reshape(filteredDDM, [], numMeas);
        dataMatrix = fillMissingWithMedian(dataMatrix);

        if numMeas == 1
            % Single-measurement fallback
            reconstructedRadargram = reshape(dataMatrix, numRows, numCols);
            reconstructedRadargram = normalizeSafe(reconstructedRadargram);
        else
            % PCA decomposition
            [~, score, ~] = pca(dataMatrix, ...
                Centered=false, VariableWeights="variance");

            principalComponent = score(:, 1);
            reconstructedRadargram = reshape(principalComponent, numRows, numCols);
        end

        reconstructedSpatialLB(:, :, cycle) = reconstructedRadargram;
        spatialMasks(:, :, cycle) = reconstructedRadargram;
    end

    %% -------------------- TEMPORAL STACK --------------------------------
    stackedDDM = reconstructedSpatialLB;

    %% -------------------- GLOBAL PCA (OPTIONAL SUMMARY) -----------------
    dataMatrix = reshape(stackedDDM, [], size(stackedDDM, 3));
    dataMatrix = fillMissingWithMedian(dataMatrix);

    [~, score, ~] = pca(dataMatrix, ...
        Centered=false, VariableWeights="variance");

    principalComponent = score(:, 1);
    reconstructedRadargram = reshape(principalComponent, numRows, numCols); %#ok<NASGU>

    %% -------------------- FFT-BASED TEMPORAL ANALYSIS -------------------
    % Annual frequency extraction from cycle stack
    N = size(stackedDDM, 3);
    dtYears = 27 / 365.25;           % Sentinel-3 repeat period in years
    Fs = 1 / dtYears;                % sampling frequency (cycles/year)

    X = stackedDDM - mean(stackedDDM, 3, 'omitnan');

    % Optional windowing
    window1D = hann(N);
    window3D = reshape(window1D, 1, 1, N);
    Xwindowed = X .* window3D;

    F = fft(Xwindowed, [], 3);

    f = (0:N-1) * (Fs / N);
    nHalf = floor(N / 2) + 1;
    fPos = f(1:nHalf);
    Fpos = F(:, :, 1:nHalf);

    % Closest frequency to 1 cycle/year
    [~, idxAnnual] = min(abs(fPos - 1));
    seasonalPattern = abs(Fpos(:, :, idxAnnual));

    %% -------------------- COMBINE SPATIAL + SEASONAL --------------------
    seasonalPatternNorm = normalizeSafe(seasonalPattern);

    for cycle = 1:numCycles
        spatialPattern = reconstructedSpatialLB(:, :, cycle);
        spatialPatternNorm = normalizeSafe(spatialPattern);

        % Combine via element-wise product
        mixedPattern = seasonalPatternNorm .* spatialPatternNorm;
        mixedPatternNorm = normalizeSafe(mixedPattern);

        % Soft mask
        spatialTemporalMasks(:, :, cycle) = mixedPatternNorm;

        % Hard thresholded version for reconstructed output
        binaryMask = mixedPattern > 0.7 * max(mixedPattern(:));
        reconstructedSpatialTemporalLB(:, :, cycle) = seasonalPattern .* binaryMask;
    end
end

%% ========================================================================
%  LOCAL FUNCTIONS
%% ========================================================================

function X = fillMissingWithMedian(X)
    validValues = X(isfinite(X));
    if isempty(validValues)
        X(:) = 0;
    else
        replacement = median(validValues, 'all');
        X(~isfinite(X)) = replacement;
    end
end

function Xnorm = normalizeSafe(X)
    maxVal = max(X(:));
    minVal = min(X(:));

    if ~isfinite(maxVal) || maxVal == 0
        Xnorm = zeros(size(X));
    elseif maxVal == minVal
        Xnorm = X ./ maxVal;
    else
        Xnorm = X ./ maxVal;
    end
end