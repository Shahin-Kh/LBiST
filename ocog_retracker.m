function [range_bin] = ocog_retracker(waveform, percentage, skip)

    waveform = waveform .^ 2;
    sq_sum = sum(waveform);
    waveform = waveform .^ 4;
    qa_sum = sum(waveform);
    
    % Compute OCOG threshold
    threshold = percentage * sqrt(qa_sum / sq_sum);
    ind_first_over = find(waveform > threshold, 1, 'first');
    
    if isempty(ind_first_over) || ind_first_over == 1
        range_bin = NaN;
        return;
    end
    
    % Linear interpolation for sub-bin accuracy
    decimal = (waveform(ind_first_over-1) - threshold) / ...
              (waveform(ind_first_over-1) - waveform(ind_first_over));
          
    range_bin = skip + ind_first_over - 1 + decimal;





