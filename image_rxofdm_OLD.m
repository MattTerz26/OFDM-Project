function [rxbits, conf] = image_rxofdm(rxsignal, conf)
    %
    %  RX for image transmission with:
    %    - SC preamble detection (frame_sync_ofdm)
    %    - periodic OFDM training symbols
    %    - Viterbi-Viterbi phase tracking
    %    - dynamic channel update
    %
    
    % ------------------ CONSTANTS --------------------------
    bitsPerSym = conf.modulation_order;      
    N          = conf.ofdm.ncarrier;         
    Ncp        = conf.ofdm.cplen;            
    fs         = conf.f_s;                   
    fc         = conf.f_c;
    os_sc      = conf.sc.os_factor;
    
    trainSym   = conf.trainSymFreq;          % training freq-domain symbol from TX
    
    % ------------------ DOWNCONVERSION ----------------------
    n = (0:length(rxsignal)-1).';
    carrier_rx = exp(-1j*2*pi*fc*n/fs);
    bb_rx = rxsignal(:).*carrier_rx;
    
    % ------------------ LOWPASS FILTER ----------------------
    bb_rx_filt = ofdmlowpass(bb_rx, conf, conf.ofdm.bandwidth);
   
    % ------------------ FRAME SYNC --------------------------
    beginning = frame_sync_ofdm_fast(bb_rx_filt, conf);
    start_ofdm_idx = beginning + 1;

    fprintf("I'm here, after frame sync");
    
    if start_ofdm_idx < 1 || start_ofdm_idx > length(bb_rx_filt)
        error('rx_task2_image: Invalid frame sync index');
    end
    
    % ------------------ REMOVE PREAMBLE ----------------------    
    ofdm_bb_os_rx = bb_rx_filt(start_ofdm_idx:end);
    
    % ------------------ RESAMPLING TO OFDM -------------------
    ofdm_bb_rx = ofdm_rx_resampling(ofdm_bb_os_rx, conf);
    
    % ------------------ SYMBOL LENGTH ------------------------
    symLen = N + Ncp;
    nSym = floor(length(ofdm_bb_rx)/symLen);
    
    if nSym < 2
        error('rx_task2_image: Not enough OFDM symbols in stream');
    end
    
    ofdm_bb_rx = ofdm_bb_rx(1 : nSym*symLen);
    rxBlocks = reshape(ofdm_bb_rx, symLen, nSym);
    
    % ------------------ REMOVE CP ----------------------------
    rxNoCP = rxBlocks(Ncp+1:end, :);
    
    % ------------------ FFT ----------------------------------
    Z = fft(rxNoCP, N, 1);    % N × nSym
    
    % ------------------ FIND TRAINING SYMBOLS -----------------
    corr_train = sum(Z .* conj(trainSym), 1);
    isTrain = abs(corr_train) > 0.7 * max(abs(corr_train));
    
    % ------------------ INITIALIZE ESTIMATION -----------------
    H = [];                          % channel estimate
    theta_hat = zeros(N,1);          % per-subcarrier phase
    alpha = 0.10;                    % smoothing factor for tracking
    
    rxbits = [];
    
    % ------------------ LOOP OVER SYMBOLS ---------------------
    for k = 1:nSym
        
        if isTrain(k)
            % -----------------CHANNEL UPDATE-----------------
            H = Z(:,k) ./ trainSym;
            theta_hat(:) = 0;      % reset phase tracker at training
            fprintf('Processed training symbol %d\n', k);
            continue;
        end
        
        if isempty(H)
            % skip symbols before first training
            continue;
        end
        
        % -------------------EQUALIZATION----------------------
        Y = Z(:,k) ./ H;
        
        % -------------------VITERBI-VITERBI PHASE TRACKING----        
        % Vectorized Viterbi-Viterbi phase tracking
        baseTheta = 0.25 * angle(-Y.^4);
        candidates = baseTheta + pi/2 * (-1:4);
        
        % Find the closest candidate for each subcarrier
        [~, idx] = min(abs(candidates - theta_hat), [], 2);
        theta_raw = candidates(idx);
        
        % Smoothing filter
        theta_hat = mod(alpha * theta_raw + (1 - alpha) * theta_hat, 2 * pi);
        
        % Correct
        Ycorr = Y .* exp(-1j * theta_hat);
        
        % ------------------- DEMAP ---------------------------
        bits_hat = demapper_QPSK(Ycorr);
        rxbits = [rxbits ; bits_hat];
        
        fprintf('Processed OFDM symbol %d\n', k);
    end

end
