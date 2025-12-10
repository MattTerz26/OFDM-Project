function [rxbits, conf] = rx_comb(rxsignal, conf)
% Receiver with:
% - 1 full training OFDM symbol
% - comb pilots in all subsequent OFDM symbols
% - LS + interpolation per symbol
% - soft phase tracking (temporal smoothing)
% - small regularization to avoid equalizer blow-up

    %% Constants
    bitsPerSym = conf.modulation_order;
    N          = conf.ofdm.ncarrier;
    Ncp        = conf.ofdm.cplen;
    fs_audio   = conf.f_s;
    fc         = conf.f_c;

    pilotIdx    = conf.ofdm.pilotIdx;
    pilotValues = conf.ofdm.pilotValues;
    dataIdx     = conf.ofdm.dataIdx;

    Ndata       = numel(dataIdx);

    nDataOfdmSym = conf.nbits / (bitsPerSym * Ndata);
    nOfdmSymTot  = nDataOfdmSym + 1;

    %% Downconversion
    n         = (0:length(rxsignal)-1).';
    carrierRx = exp(-1j * 2*pi*fc * n / fs_audio);
    bb_rx     = rxsignal(:) .* carrierRx;

    %% Low-pass filtering
    bb_rx_filt = ofdmlowpass(bb_rx, conf, conf.ofdm.bandwidth);

    %% Frame sync
    beginning_of_data = frame_sync_ofdm(bb_rx_filt, conf);
    start_ofdm_idx    = beginning_of_data + 1;

    if start_ofdm_idx <= 0 || start_ofdm_idx > length(bb_rx_filt)
        error('rx_comb: start index out of range (%d).', start_ofdm_idx);
    end

    ofdm_bb_os_rx = bb_rx_filt(start_ofdm_idx:end);

    %% Resample to OFDM Fs
    ofdm_bb_rx = ofdm_rx_resampling(ofdm_bb_os_rx, conf);

    %% Trim to expected OFDM length
    totalLen = nOfdmSymTot * (N + Ncp);
    if length(ofdm_bb_rx) < totalLen
        error('rx_comb: signal too short (%d < %d).', length(ofdm_bb_rx), totalLen);
    end
    ofdm_bb_rx = ofdm_bb_rx(1:totalLen);

    %% Reshape + remove CP
    rxBlocks = reshape(ofdm_bb_rx, N+Ncp, nOfdmSymTot);
    rxNoCP   = rxBlocks(Ncp+1:end, :);

    %% FFT
    Z = fft(rxNoCP, N, 1);     % columns: [training | data1 | data2 ...]

    %% Full training (symbol 1)
    Z_train      = Z(:,1);
    trainSymFreq = conf.ofdm.trainSymFreq;
    H0           = Z_train ./ trainSymFreq;

    %% Data symbols
    Z_data   = Z(:,2:end);
    nDataSym = size(Z_data, 2);

    H_est  = zeros(N, nDataSym);
    Y      = zeros(N, nDataSym);
    Y_corr = zeros(N, nDataSym);

    k_all = (1:N).';

    %% Smoothing parameters
    alphaH   = 0.3;     % temporal smoothing of pilot-based LS estimate
    alphaPhi = 0.3;     % smoothing of residual phase
    epsH     = 1e-3;    % avoid |H| → 0 in deep fades

    %% Initial values from training symbol
    H_p_prev = H0(pilotIdx);
    phi_prev = 0;

    %% Process each OFDM data symbol
    for s = 1:nDataSym
        
        Zs = Z_data(:, s);

        %% Pilot-based LS estimation
        Z_p_raw = Zs(pilotIdx);
        H_p_raw = Z_p_raw ./ pilotValues;

        %% Temporal smoothing of pilot channel estimate
        H_p = alphaH * H_p_raw + (1 - alphaH) * H_p_prev;
        H_p_prev = H_p;

        %% Frequency interpolation (shape-preserving cubic)
        H_est(:, s) = interp1(pilotIdx(:), H_p(:), k_all, 'pchip', 'extrap');

        %% Equalization with a small regularization term
        Hs = H_est(:, s);
        Hs_mag = abs(Hs);

        small = (Hs_mag < epsH);
        if any(small)
            Hs(small) = epsH .* exp(1j * angle(Hs(small)));
        end

        Y(:, s) = Zs ./ Hs;

        %% Soft phase tracking (global per symbol)
        Y_pilot_eq = Y(pilotIdx, s) ./ pilotValues;
        phi_meas   = angle(sum(Y_pilot_eq));

        phi_soft = angle(exp(1j * (alphaPhi * phi_meas + (1 - alphaPhi) * phi_prev)));
        phi_prev = phi_soft;

        Y_corr(:, s) = Y(:, s) * exp(-1j * phi_soft);
    end

    %% Extract only data carriers
    A_data = Y_corr(dataIdx, :);
    A_vec  = A_data(:);

    rxbits = demapper_QPSK(A_vec);

    if length(rxbits) < conf.nbits
        error('rx_comb: not enough output bits (%d < %d).', length(rxbits), conf.nbits);
    end

    rxbits = rxbits(1:conf.nbits);
end
