% Number of OFDM data symbols (excluding training)
    nDataOfdmSym = conf.nbits / (bitsPerSym * N);
    nOfdmSymTot  = nDataOfdmSym + 1;           % +1 training symbol

    % Remove leading/trailing zeros
    rx_trim = rxsignal(fs_audio+1 : end-fs_audio);   % keep only txsignal part

    % Downconversion to complex baseband
    n = (0:length(rx_trim)-1).';
    carrier_rx = exp(-1j*2*pi*fc*n/fs_audio);       % minus sign for RX
    bb_rx = rx_trim(:) .* carrier_rx;

    % Low-pass filtering (extract baseband)
    bb_rx_filt = ofdmlowpass(bb_rx, conf, conf.ofdm.bandwidth);

    % Remove preamble in time-domain by knowing its length
    %    Preamble length:
    %       - conf.sc.nsyms symbols
    %       - upsampled by os_sc
    %       - convolved with RRC of length Lp = 2*txpulse_length + 1
    %       => Np = nsyms*os_sc + Lp - 1

    Lp = length(conf.sc.txpulse);
    preamble_len = conf.sc.nsyms*os_sc + Lp - 1;

    if length(bb_rx_filt) <= preamble_len
        error('rx_task1_basic: received signal too short compared to preamble length.');
    end

    ofdm_bb_os_rx = bb_rx_filt(preamble_len+1:end);

    % Resample from audio Fs to OFDM Fs
    ofdm_bb_rx = ofdm_rx_resampling(ofdm_bb_os_rx, conf);

    % Shorten to expected length (multiple of (N+Ncp)*nOfdmSymTot )
    total_ofdm_len = nOfdmSymTot * (N + Ncp);

    if length(ofdm_bb_rx) < total_ofdm_len
        error('rx_task1_basic: ofdm_bb_rx too short (got %d, need %d).', ...
               length(ofdm_bb_rx), total_ofdm_len);
    end
    ofdm_bb_rx = ofdm_bb_rx(1:total_ofdm_len);