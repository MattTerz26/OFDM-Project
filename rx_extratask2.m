function [rxbits, conf, H_est, Z_data, Y_corr] = rx_extratask2(rxsignal, conf)
% RX_TASK2 - Receiver for Audio OFDM with Scrambling and Phase Tracking

    % Constants
    bitsPerSym   = conf.modulation_order;       
    N            = conf.ofdm.ncarrier;          
    Ncp          = conf.ofdm.cplen;             
    fs_audio     = conf.f_s;                    
    fc           = conf.f_c;                    
    
    % Derived lengths
    nDataOfdmSym = conf.nbits / (bitsPerSym * N);
    nOfdmSymTot  = nDataOfdmSym + 1;           % +1 training symbol
    
    % 1. Downconversion
    n = (0:length(rxsignal)-1).';
    carrier_rx = exp(-1j*2*pi*fc*n/fs_audio);
    bb_rx = rxsignal(:) .* carrier_rx;
    
    % 2. Low-pass filtering
    bb_rx_filt = ofdmlowpass(bb_rx, conf, conf.ofdm.bandwidth);
    
    % 3. Synchronization (Using YOUR correct function)
    beginning_of_data = frame_sync_ofdm_fast(bb_rx_filt, conf);        
    
    start_ofdm_idx = beginning_of_data + 1;
    
    % Safety check
    if start_ofdm_idx <= 0 || start_ofdm_idx > length(bb_rx_filt)
        warning('rx_task2: Sync failed or index out of bounds. Defaulting to 1.');
        start_ofdm_idx = 1;
    end
    
    % 4. Extraction
    ofdm_bb_os_rx = bb_rx_filt(start_ofdm_idx:end);
    
    % 5. Resample
    ofdm_bb_rx = ofdm_rx_resampling(ofdm_bb_os_rx, conf);
    
    % 6. Cut to exact length
    total_ofdm_len = nOfdmSymTot * (N + Ncp);
    if length(ofdm_bb_rx) < total_ofdm_len
        % Zero-pad if recording cut off too early
        padding = zeros(total_ofdm_len - length(ofdm_bb_rx), 1);
        ofdm_bb_rx = [ofdm_bb_rx; padding];
    end
    ofdm_bb_rx = ofdm_bb_rx(1:total_ofdm_len);
    
    % 7. Parallelize & Remove CP
    rxBlocks = reshape(ofdm_bb_rx, N+Ncp, nOfdmSymTot);
    rxNoCP = rxBlocks(Ncp+1:end, :);    
    
    % 8. FFT
    Z = fft(rxNoCP, N, 1);              
    
    % 9. Channel Estimation (Training Symbol)
    Z_train = Z(:,1);                   
    trainSymFreq = conf.ofdm.trainSymFreq;   
    H_est = Z_train ./ trainSymFreq;
    
    % 10. Equalization
    Z_data = Z(:,2:end);               
    nDataSym = size(Z_data,2);
    
    H_mat = repmat(H_est, 1, nDataSym);
    Y = Z_data ./ H_mat;          
    
    % 11. Phase Tracking (FIXED)
    if ~isfield(conf.ofdm, 'alpha')
        alpha = 0.1;
    else
        alpha = conf.ofdm.alpha;
    end
    
    theta_hat = zeros(N, nDataSym+1);
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
    
    % 12. Demapping
    A_hat = Y_corr(:);
    rxbits = demapper_QPSK(A_hat);
    
    % Trim to expected bits
    if length(rxbits) > conf.nbits
        rxbits = rxbits(1:conf.nbits);
    end
    
    % 13. Descrambling (FIXED Variable Name)
    % You cannot use 'length(bits_in)' here because bits_in doesn't exist.
    % We use 'conf.nbits' or the current length.
    
    if isfield(conf, 'enable_scrambler') && conf.enable_scrambler
        rng(42); % Must match TX seed
        scramble_seq = randi([0 1], length(rxbits), 1);
        rxbits = xor(rxbits, scramble_seq);
    end
end