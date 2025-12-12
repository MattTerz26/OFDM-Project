function [rxbits, conf] = rxofdm(rxsignal, conf)
% Digital Receiver (dispatcher)
%
%   [rxbits, conf] = rxofdm(rxsignal, conf) implements different receiver
%   variants depending on conf.rx_mode. It serves as a central function
%   that routes the received signal to the appropriate processing function
%   based on the specified receiver mode in the configuration structure.
%
%   Inputs:
%       rxsignal : received real passband signal (from emulator/audio)
%       conf     : configuration structure containing receiver settings
%
%   Outputs:
%       rxbits   : recovered information bits from the received signal
%       conf     : updated configuration structure (if modified)
%
%   Receiver modes:
%       - "task1": Calls the rx_task1 function for the simplest OFDM receiver
%       - "task2": Calls the rx_task2 function for the OFDM receiver with phase tracking
%       - "extratask1": Calls the rx_extratask1 function for images /
%       packets transmission
%
%   Error Handling:
%       If an unknown receiver mode is specified, an error is raised indicating
%       the unknown mode.

    % Default receiver mode if not specified
    if ~isfield(conf, 'rx_mode')
        conf.rx_mode = "task1";
    end

    switch string(conf.rx_mode)
        case "task1"
            [rxbits, conf] = rx_task1(rxsignal, conf);

        case "task2"
            [rxbits, conf] = rx_task2(rxsignal, conf);

        case "extratask1"
            [rxbits, conf] = rx_extratask1(rxsignal, conf);

        case "combtype"
            [rxbits, conf] = rx_comb_alt(rxsignal, conf);

        otherwise
            error("rxofdm: Unknown receiver mode '%s'", conf.rx_mode);
    end
end

%% ==================================================================
%  Receiver for Task 1: simplest OFDM receiver (no phase tracking)
%  ==================================================================
function [rxbits, conf] = rx_task1(rxsignal, conf)

    % Useful constants
    bitsPerSym   = conf.modulation_order;       % 2 for QPSK
    N            = conf.ofdm.ncarrier;          % 512
    Ncp          = conf.ofdm.cplen;             % 256
    fs_audio     = conf.f_s;                    % 48 kHz
    fc           = conf.f_c;                    % 8 kHz
    os_sc        = conf.sc.os_factor;           % oversampling preamble

    % Number of OFDM data symbols (excluding training)
    nDataOfdmSym = conf.nbits / (bitsPerSym * N);
    nOfdmSymTot  = nDataOfdmSym + 1;           % +1 training symbol

    % Downconversion to complex baseband    - Time domain
    n = (0:length(rxsignal)-1).';
    carrier_rx = exp(-1j*2*pi*fc*n/fs_audio);       % minus sign for RX
    bb_rx = rxsignal(:) .* carrier_rx;

    % Low-pass filtering (extract baseband) - Time domain
    bb_rx_filt = ofdmlowpass(bb_rx, conf, conf.ofdm.bandwidth);

    % filtering for preamble
    rx_filt_p_det = lowpass(real(bb_rx), 2000, conf.f_s);
    
    beginning_of_data = frame_sync_ofdm(rx_filt_p_det, conf);

    start_ofdm_idx = beginning_of_data + 1;

    if start_ofdm_idx <= 0 || start_ofdm_idx > length(bb_rx_filt)
        error('rx_task1: start_ofdm_idx out of range (%d).', start_ofdm_idx);
    end

    % OFDM signal
    ofdm_bb_os_rx = bb_rx_filt(start_ofdm_idx:end);

    % Resample from audio Fs to OFDM Fs
    ofdm_bb_rx = ofdm_rx_resampling(ofdm_bb_os_rx, conf);

    % Shorten to expected length (multiple of (N+Ncp)*nOfdmSymTot )
    total_ofdm_len = nOfdmSymTot * (N + Ncp);

    if length(ofdm_bb_rx) < total_ofdm_len
        error('rx_task1: ofdm_bb_rx too short (got %d, need %d).', ...
               length(ofdm_bb_rx), total_ofdm_len);
    end
    ofdm_bb_rx = ofdm_bb_rx(1:total_ofdm_len);

    % Parallelize
    rxBlocks = reshape(ofdm_bb_rx, N+Ncp, nOfdmSymTot);

    % Remove CP
    rxNoCP = rxBlocks(Ncp+1:end, :);    % size N x nOfdmSymTot

    % FFT to frequency domain (vector form)
    Z = fft(rxNoCP, N, 1);              % columns: [training | data ...]

    % Channel estimation from training symbol
    Z_train = Z(:,1);                   % training symbol in freq
    trainSymFreq = conf.ofdm.trainSymFreq;   % saved in TX
    
    % Channel estimation
    H_est = Z_train ./ trainSymFreq;

    % Equalization and demapping of data symbols
    Z_data = Z(:,2:end);               % all data OFDM symbols
    nDataSym = size(Z_data,2);

    % Equalize
    H_mat  = repmat(H_est, 1, nDataSym);
    A_hat  = Z_data ./ H_mat;          % estimated QPSK symbols in freq

    % Demap
    rxbits = demapper_QPSK(A_hat);
    
    % Keep only the first conf.nbits bits (safety)
    if length(rxbits) < conf.nbits
        error('rx_task1: not enough bits decoded (got %d, need %d).', ...
               length(rxbits_all), conf.nbits);
    end
    rxbits = rxbits(1:conf.nbits);

end

%% ==================================================================
%  Receiver for Task 2: phase tracking
% ===================================================================
function [rxbits, conf]  = rx_task2(rxsignal, conf)

    % Useful constants
    bitsPerSym   = conf.modulation_order;       % 2 for QPSK
    N            = conf.ofdm.ncarrier;          % 512
    Ncp          = conf.ofdm.cplen;             % 256
    fs_audio     = conf.f_s;                    % 48 kHz
    fc           = conf.f_c;                    % 8 kHz
    os_sc        = conf.sc.os_factor;           % oversampling preamble

    % Number of OFDM data symbols (excluding training)
    nDataOfdmSym = conf.nbits / (bitsPerSym * N);
    nOfdmSymTot  = nDataOfdmSym + 1;           % +1 training symbol

    % Downconversion to complex baseband    - Time domain
    n = (0:length(rxsignal)-1).';
    carrier_rx = exp(-1j*2*pi*fc*n/fs_audio);       % minus sign for RX
    bb_rx = rxsignal(:) .* carrier_rx;

    % Low-pass filtering (extract baseband) - Time domain
    bb_rx_filt = ofdmlowpass(bb_rx, conf, conf.ofdm.bandwidth);
    
    beginning_of_data = frame_sync_ofdm_fast(bb_rx_filt, conf);        

    start_ofdm_idx = beginning_of_data + 1;

    if start_ofdm_idx <= 0 || start_ofdm_idx > length(bb_rx_filt)
        error('rx_task2: start_ofdm_idx out of range (%d).', start_ofdm_idx);
    end

    % OFDM signal
    ofdm_bb_os_rx = bb_rx_filt(start_ofdm_idx:end);

    % Resample from audio Fs to OFDM Fs
    ofdm_bb_rx = ofdm_rx_resampling(ofdm_bb_os_rx, conf);

    % Shorten to expected length (multiple of (N+Ncp)*nOfdmSymTot )
    total_ofdm_len = nOfdmSymTot * (N + Ncp);

    if length(ofdm_bb_rx) < total_ofdm_len
        error('rx_task1: ofdm_bb_rx too short (got %d, need %d).', ...
               length(ofdm_bb_rx), total_ofdm_len);
    end
    ofdm_bb_rx = ofdm_bb_rx(1:total_ofdm_len);

    % Parallelize
    rxBlocks = reshape(ofdm_bb_rx, N+Ncp, nOfdmSymTot);

    % Remove CP
    rxNoCP = rxBlocks(Ncp+1:end, :);    % size N x nOfdmSymTot

    % FFT to frequency domain (vector form)
    Z = fft(rxNoCP, N, 1);              % columns: [training | data ...]

    % Channel estimation from training symbol
    Z_train = Z(:,1);                   % training symbol in freq
    trainSymFreq = conf.ofdm.trainSymFreq;   % saved in TX

    % Channel estimation
    H_est = Z_train ./ trainSymFreq;

    % Equalization and demapping of data symbols
    Z_data = Z(:,2:end);               % all data OFDM symbols
    nDataSym = size(Z_data,2);

    % Equalize
    H_mat  = repmat(H_est, 1, nDataSym);
    Y  = Z_data ./ H_mat;          % estimated QPSK symbols in freq
    
    % Phase tracking and phase correction
    
    if ~isfield(conf.ofdm, 'alpha')
        alpha = 0.1;  % Default value if alpha is not defined
    else
        alpha = conf.ofdm.alpha;  % IIR filter coefficient
    end

    theta_hat = zeros(N, nDataSym+1);
    theta_hat(:, 1) = 0;

    Y_corr = zeros(size(Y));
    
    for n = 1:nDataSym      % looping through OFDM data symbols
        y = Y(:,n);         % vectorized access to all carriers for the current symbol
        
        % Viterbi-Viterbi
        baseTheta = 0.25 * angle(-y.^4);
        candidates = baseTheta + pi/2 * (-1:4);   % 6 candidates
        
        % closest candidate
        [~, idx] = min(abs(candidates - theta_hat(:, n)), [], 2);
        theta_raw = candidates(sub2ind(size(candidates), (1:N).', idx(:))); % vectorized indexing

        % IIR filter
        theta_hat(:, n+1) = alpha*theta_raw + (1-alpha)*theta_hat(:, n);
        theta_hat(:, n+1) = mod(theta_hat(:, n+1), 2*pi);

        % final correction
        Y_corr(:,n) = y .* exp(-1j*theta_hat(:, n+1));

    end
    
    A_hat = Y_corr(:);

    % Demap
    rxbits = demapper_QPSK(A_hat);
    
    % Keep only the first conf.nbits bits (safety)
    if length(rxbits) < conf.nbits
        error('rx_task2: not enough bits decoded (got %d, need %d).', ...
               length(rxbits_all), conf.nbits);
    end
    rxbits = rxbits(1:conf.nbits);

end

%% ==================================================================
%  Receiver for Extra Task 1: image transmission
% ===================================================================
function [rxbits, conf] = rx_extratask1(rxsignal, conf)
    %
    %  RX for image transmission with:
    %    - SC preamble detection
    %    - periodic OFDM training symbols
    %    - Viterbi-Viterbi phase tracking
    %    - dynamic channel update
    %
    
    % Constants
    bitsPerSym = conf.modulation_order;      
    N          = conf.ofdm.ncarrier;         
    Ncp        = conf.ofdm.cplen;            
    fs         = conf.f_s;                   
    fc         = conf.f_c;
    os_sc      = conf.sc.os_factor;
    
    trainSym   = conf.trainSymFreq;          % Training frequency-domain symbol from TX
    
    % Downconversion
    n = (0:length(rxsignal)-1).';
    carrier_rx = exp(-1j*2*pi*fc*n/fs);
    bb_rx = rxsignal(:).*carrier_rx;
    
    % Low-pass filtering
    bb_rx_filt = ofdmlowpass(bb_rx, conf, conf.ofdm.bandwidth);
   
    % filtering for preamble
    rx_filt_p_det = lowpass(real(bb_rx_filt), 2000, conf.f_s);

    % Frame synchronization
    beginning = frame_sync_ofdm_fast(rx_filt_p_det, conf);
    start_ofdm_idx = beginning + 1;

    if start_ofdm_idx < 1 || start_ofdm_idx > length(bb_rx_filt)
        error('ExtraTask1: Invalid frame sync index');
    end
    
    % Remove preamble    
    ofdm_bb_os_rx = bb_rx_filt(start_ofdm_idx:end);
    
    % Resampling to OFDM
    ofdm_bb_rx = ofdm_rx_resampling(ofdm_bb_os_rx, conf);
    
    % Symbol length
    symLen = N + Ncp;
    nSym = floor(length(ofdm_bb_rx)/symLen);
    
    if nSym < 2
        error('ExtraTask1: Not enough OFDM symbols in stream');
    end
    
    ofdm_bb_rx = ofdm_bb_rx(1 : nSym*symLen);

    % Parallelize
    rxBlocks = reshape(ofdm_bb_rx, symLen, nSym);
    
    % Remove CP
    rxNoCP = rxBlocks(Ncp+1:end, :);
    
    % FFT
    Z = fft(rxNoCP, N, 1);    % N × nSym
    
    % Find training symbols
    corr_train = sum(Z .* conj(trainSym), 1);
    isTrain = abs(corr_train) > 0.7 * max(abs(corr_train));
    
    % Initialize estimation
    % --- INITIALIZATION ---
    H_est_amp = [];                 % amplitude-only channel estimate
    theta_hat = zeros(N,1);         % per-carrier tracked phase
    haveH     = false;
    
    alpha = 0.1;                   % smoothing for data symbols
    
    rxbits = [];
    
    % --- SYMBOL LOOP ---
    for k = 1:nSym
        
        % ===== TRAINING SYMBOL =====
        if isTrain(k)
            % Estimate channel (complex)
            H_temp = Z(:,k) ./ trainSym;
    
            % Amplitude-only channel estimate
            H_est_amp = abs(H_temp);
    
            % Initialize per-carrier phase tracker
            theta_hat = angle(H_temp);    % N×1
    
            fprintf("Processed training symbol %d\n",k);
            continue;
        end
    
        % ===== DATA SYMBOL =====
    
        % 1) Amplitude equalization ONLY
        Y_amp = Z(:,k) ./ H_est_amp;      % N×1
    
        % 2) Compute Viterbi–Viterbi raw phase (per carrier)
        delta = 0.25 * angle(-Y_amp.^4);  % N×1
    
        % 3) Build rotation candidates (N×6)
        rotations = pi/2 * (-1:4);        % 1×6
        candidates = delta + rotations;   % N×6
    
        % 4) Compare each candidate with previous theta_hat
        prev_rep = theta_hat .* ones(1,6);     % N×6
        [~, idx] = min(abs(candidates - prev_rep), [], 2);
    
        % 5) Extract chosen candidate for each carrier
        theta_new = candidates(sub2ind(size(candidates),(1:N).', idx));
    
        % 6) Smooth per-carrier phase (IIR tracking)
        theta_hat = alpha * theta_new + (1-alpha) * theta_hat;
    
        % 7) Apply phase correction (AFTER amp equalization)
        Y_corr = Y_amp .* exp(-1j * theta_hat);
    
        % 8) Demap
        bits_hat = demapper_QPSK(Y_corr(:));
        rxbits   = [rxbits ; bits_hat];

    end
    
    % Keep only the first conf.nbits bits (safety)
    if length(rxbits) < conf.nbits
        error('rx_task2: not enough bits decoded (got %d, need %d).', ...
               length(rxbits_all), conf.nbits);
    end
    rxbits = rxbits(1:conf.nbits);

    % Check for scrambling configuration
    if isfield(conf.ofdm, 'scramble_seq')
        scramble_seq = conf.ofdm.scramble_seq;
    else
        rng(42);
        scramble_seq = randi([0 1], length(bits_in), 1);
    end
    rxbits = xor(rxbits, scramble_seq);

end

%% ==================================================================
%  Receiver for comb type training
% ===================================================================

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

    bitsPerOfdmDataSym = bitsPerSym * Ndata;
    nDataOfdmSym = conf.nbits_padded / bitsPerOfdmDataSym;
    if mod(conf.nbits_padded, bitsPerOfdmDataSym) ~= 0
        error('rx_comb: padded bit count not multiple of OFDM payload');
    end

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
