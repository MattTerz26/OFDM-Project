function beginning_of_data = frame_sync_ofdm_fast(rx_bb, conf)
    % FRAME_SYNC_OFDM_FAST - Synchronizes the received OFDM signal
    % 
    % Syntax: beginning_of_data = frame_sync_ofdm_fast(rx_bb, conf)
    %
    % Inputs:
    %   rx_bb - Received baseband signal (vector)
    %   conf - Configuration structure containing:
    %       conf.f_s - Sampling frequency
    %       conf.sc.preamble_bb - Expected preamble signal (vector)
    %
    % Outputs:
    %   beginning_of_data - Index of the beginning of the data in the received signal
    %
    % Description:
    % This function performs frame synchronization for an OFDM signal by
    % correlating the received signal with the expected preamble. It uses
    % fast correlation through filtering to identify the start of the data.
    
    % Expected preamble region
    margin = 3000;
    expected = conf.f_s;
    search_end   = min(length(rx_bb), expected + margin + length(conf.sc.preamble_bb));

    rx_sub = rx_bb(1:end);

    % Preamble
    p = conf.sc.preamble_bb(:);
    Lp = length(p);

    % Fast correlation using filter
    c = filter(conj(flipud(p)), 1, rx_sub);
    pow_rx = abs(rx_sub).^2;
    e = filter(ones(length(p),1), 1, pow_rx);

    T = abs(c).^2 ./ (e + eps);                       % normalized statistic

    % ignore the initial part where the window is not full
    [~, idx] = max(T(Lp:end));
    beginning_of_data = idx + Lp - 1;    
end