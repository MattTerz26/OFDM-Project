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
            rxbits = rx_task1(rxsignal, conf);

        case "task2"
            rxbits = rx_task2(rxsignal, conf);

        case "extratask1"
            rxbits = rx_extratask1(rxsignal, conf);

        otherwise
            error("rxofdm: Unknown receiver mode '%s'", conf.rx_mode);
    end
end

%% ==================================================================
%  Receiver for Task 1: simplest OFDM receiver (no phase tracking)
%  ==================================================================
function rxbits = rx_task1(rxsignal, conf)

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
function rxbits = rx_task2(rxsignal, conf)

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
    theta_hat = zeros(N, nDataSym+1);
    alpha = 1;                 % IIR filter coefficient
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
        theta_hat(:, n+1) = mod(alpha*theta_raw + (1-alpha)*theta_hat(:, n), 2*pi);
        
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
    H_est = [];                          % Channel estimate
    theta_hat = zeros(N,nSym+1);          % Per-subcarrier phase
    alpha = 0.01;                    % Smoothing factor for tracking
    
    rxbits = [];
    
    % Loop over symbols
    for k = 1:nSym
        
        if isTrain(k)
            % Channel update
            H_est = Z(:,k) ./ trainSym;
            theta_hat(:,k) = 0;      % Reset phase tracker at training
            % fprintf('Processed training symbol %d\n', k);
            continue;
        end
        
        if isempty(H_est)
            % Skip symbols before first training
            continue;
        end
        
        % Equalization
        Y = Z(:,k) ./ H_est;
        
        % Viterbi-Viterbi phase tracking        
        baseTheta = 0.25 * angle(-Y.^4);
        candidates = baseTheta + pi/2 * (-1:4);
        
        % Find the closest candidate for each subcarrier
        [~, idx] = min(abs(candidates - theta_hat(:, k-1)), [], 2);
        theta_raw = candidates(idx);
        
        % IIR filter
        theta_hat(:, k) = mod(alpha * theta_raw + (1 - alpha) * theta_hat(:, k-1), 2 * pi);
        
        % Correction
        Y_corr = Y .* exp(-1j * theta_hat(:, k));
        
        % Demapping
        bits_hat = demapper_QPSK(Y_corr(:));
        rxbits = [rxbits ; bits_hat];
        
        fprintf('Processed OFDM symbol %d\n', k);
    end

end
