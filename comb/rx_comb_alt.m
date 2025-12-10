function [rxbits, conf] = rx_comb_alt(rxsignal, conf)
% Receiver with:
% - 1 full training OFDM symbol
% - diagonal moving comb pilots (checkerboard pattern)
% - LS + interpolation, temporal smoothing, soft phase tracking

    bitsPerSym = conf.modulation_order;
    N          = conf.ofdm.ncarrier;
    Ncp        = conf.ofdm.cplen;
    fs_audio   = conf.f_s;
    fc         = conf.f_c;

    pilotSpacing = conf.ofdm.pilotSpacing;
    baseIdx      = conf.ofdm.pilotBaseIdx;
    pilotShift   = conf.ofdm.pilotShift;
    pilotValues  = conf.ofdm.pilotValues;

    Npil  = numel(baseIdx);
    Ndata = N - Npil;

    nDataOfdmSym = conf.nbits / (bitsPerSym * Ndata);
    nOfdmSymTot  = nDataOfdmSym + 1;

    %% Downconversion
    n = (0:length(rxsignal)-1).';
    carrier_rx = exp(-1j * 2*pi*fc*n/fs_audio);
    bb_rx      = rxsignal(:) .* carrier_rx;

    %% Low-pass filtering
    bb_rx_filt = ofdmlowpass(bb_rx, conf, conf.ofdm.bandwidth);

    %% Frame sync
    beginning_of_data = frame_sync_ofdm(bb_rx_filt, conf);
    start_ofdm_idx    = beginning_of_data + 1;

    if start_ofdm_idx <= 0 || start_ofdm_idx > length(bb_rx_filt)
        error('rx_comb: start index out of range (%d).', start_ofdm_idx);
    end

    ofdm_bb_os_rx = bb_rx_filt(start_ofdm_idx:end);

    %% Resample to OFDM sampling rate
    ofdm_bb_rx = ofdm_rx_resampling(ofdm_bb_os_rx, conf);

    %% Trim to expected length
    total_ofdm_len = nOfdmSymTot * (N + Ncp);
    if length(ofdm_bb_rx) < total_ofdm_len
        error('rx_comb: signal too short (%d, need %d).', ...
              length(ofdm_bb_rx), total_ofdm_len);
    end
    ofdm_bb_rx = ofdm_bb_rx(1:total_ofdm_len);

    %% Reshape + remove CP
    rxBlocks = reshape(ofdm_bb_rx, N+Ncp, nOfdmSymTot);
    rxNoCP   = rxBlocks(Ncp+1:end, :);

    %% FFT
    Z = fft(rxNoCP, N, 1);

    %% Full training symbol
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
    alphaH   = 0.3;
    alphaPhi = 0.3;
    epsH     = 1e-3;

    %% Initial state
    H_p_prev = H0(baseIdx);
    phi_prev = 0;

    A_vec = zeros(Ndata * nDataSym, 1);
    ptr   = 1;

    %% Process each OFDM data symbol
    for s = 1:nDataSym

        Zs = Z_data(:, s);

        % diagonal pilot index for symbol s
        shift        = (s-1)*pilotShift;
        pilotIdx_s   = mod(baseIdx - 1 + shift, N) + 1;
        pilotIdx_s   = sort(pilotIdx_s);

        dataIdx_s = setdiff(1:N, pilotIdx_s);
        if numel(dataIdx_s) ~= Ndata
            error('rx_comb: wrong dataIdx_s size (%d != %d) in sym %d.', ...
                  numel(dataIdx_s), Ndata, s);
        end

        %% LS channel estimation on pilots
        Z_p_raw = Zs(pilotIdx_s);
        H_p_raw = Z_p_raw ./ pilotValues;

        %% Temporal smoothing
        H_p = alphaH * H_p_raw + (1 - alphaH) * H_p_prev;
        H_p_prev = H_p;

        %% Frequency interpolation
        H_est(:, s) = interp1(pilotIdx_s(:), H_p(:), k_all, 'pchip', 'extrap');

        %% Equalization (with regularization)
        Hs = H_est(:, s);
        Hs_mag = abs(Hs);
        small = (Hs_mag < epsH);
        if any(small)
            Hs(small) = epsH .* exp(1j * angle(Hs(small)));
        end

        Y(:, s) = Zs ./ Hs;

        %% Soft phase tracking (one global phase per symbol)
        Y_pilot_eq = Y(pilotIdx_s, s) ./ pilotValues;
        phi_meas   = angle(sum(Y_pilot_eq));

        phi_soft = angle(exp(1j * (alphaPhi * phi_meas + (1 - alphaPhi) * phi_prev)));
        phi_prev = phi_soft;

        Y_corr(:, s) = Y(:, s) * exp(-1j * phi_soft);

        %% Collect data subcarriers
        dataSym_s = Y_corr(dataIdx_s, s);
        A_vec(ptr : ptr+Ndata-1) = dataSym_s;
        ptr = ptr + Ndata;
    end

    %% Demap to bits
    rxbits = demapper_QPSK(A_vec);

    if length(rxbits) < conf.nbits
        error('rx_comb: not enough bits (%d < %d).', length(rxbits), conf.nbits);
    end

    rxbits = rxbits(1:conf.nbits);
end
