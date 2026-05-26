function retrack_bin = Retracker4_Thereshold(W, threshold , l)
    % THRESHOLD_RETRACKING - Implements the Threshold Retracking Algorithm.
    % Inputs:
    %   waveform  - 1D array of waveform power values
    %   threshold - Threshold percentage (e.g., 20 for 20%)
    % Output:
    %   retrack_bin - Interpolated range bin index of the retracking point

    % Ensure waveform is a column vector
    % waveform = waveform(:);
    
    % Find the maximum waveform power
    % max_power = max(waveform);
    [COG,Amp,Width,rtrckr_bin]=Retracker1_OCOG(W,20,l);
    % Compute the threshold power level
    threshold_value = (threshold / 100) * Amp;
    
    % Find the first bin where waveform exceeds the threshold
    idx1 = find(W >= threshold_value, 1, 'first');

    % If no valid point is found, return NaN
    if isempty(idx1) || idx1 == 1
        warning('No valid retracking point found. Check waveform or threshold value.');
        retrack_bin = NaN;
        return;
    end

    % Perform linear interpolation for better accuracy
    idx0 = idx1 - 1;
    x0 = idx0; x1 = idx1;
    y0 = W(idx0); y1 = W(idx1);
    retrack_bin = x0 + (threshold_value - y0) * (x1 - x0) / (y1 - y0);
end
