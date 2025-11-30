function beginning_of_data = frame_sync_ofdm(rx_bb, conf)
%FRAME_SYNC_OFDM
%  rx_bb                : complex baseband signal (after downconversion + low-pass)
%  conf.sc.preamble_bb  : COMPLEX waveform of the preamble (as in TX)
%  beginning_of_data    : index (1-based) of the LAST sample of the preamble
%                        in the rx_bb vector

    p = conf.sc.preamble_bb(:);          % known preamble waveform (48 kHz)
    Lp = length(p);

    if length(rx_bb) < Lp
        error('frame_sync_ofdm: rx_bb too short (%d samples, preamble = %d)', ...
              length(rx_bb), Lp);
    end

    % Parameters for your frame_sync style
    detection_threshold   = 15;    % you can tune if needed
    current_peak_value    = 0;
    samples_after_threshold = Lp;  % search window after first threshold crossing

    found = false;
    beginning_of_data = -1;

    % Loop over time and perform normalized "sliding" correlation
    for i = Lp : length(rx_bb)
        seg = rx_bb(i-Lp+1 : i);        % window of Lp samples
        c   = p' * seg;                 % complex correlation
        T   = abs(c)^2 / (seg' * seg);  % normalized statistic

        if (T > detection_threshold || samples_after_threshold < Lp)
            samples_after_threshold = samples_after_threshold - 1;

            if (T > current_peak_value)
                current_peak_value  = T;
                beginning_of_data   = i;        % last sample of the preamble
            end

            if samples_after_threshold == 0
                found = true;
                break;
            end
        end
    end

    if ~found || beginning_of_data < 0
        error('frame_sync_ofdm: preamble not found (threshold too high?)');
    end
end
