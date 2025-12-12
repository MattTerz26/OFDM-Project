function [txsignal, conf] = image_txofdm(bits_in, conf)
    % IMAGE_TXOFDM - Transmitter for image packets with periodic training OFDM
    %
    % bits_in : full bitstream of the image
    % conf    : system configuration
    %
    % Output:
    %   txsignal : real-valued audio waveform ready for playback
    %
    
    % Parameters
    N      = conf.ofdm.ncarrier;
    Ncp    = conf.ofdm.cplen; 
    bps    = conf.modulation_order;      % 2 for QPSK
    bitsPerOFDM = N * bps;
    fs_audio    = conf.f_s;                  
    fc          = conf.f_c;                  
    
    % Training periodicity
    if ~isfield(conf, "train_period")
        conf.train_period = 20;   % every 20 data-OFDM symbols
    end
    
    % Packetization
    nbits = length(bits_in);
    nDataOFDM = ceil(nbits / bitsPerOFDM);
    padding = nDataOFDM * bitsPerOFDM - nbits;
    
    bits_tx = [bits_in ; zeros(padding, 1)];
    conf.nbits_padded = length(bits_tx);
    
    % Reshape into packets
    data_bits_mat = reshape(bits_tx, bitsPerOFDM, nDataOFDM);
    
    % Training OFDM
    % Generate random training bits
    trainBits = preamble_generate(bitsPerOFDM);
    trainSym  = mapper_QPSK(trainBits);
    
    conf.trainSymFreq = trainSym;
    
    % Build OFDM Sequence
    % Convert data bits to QPSK for all data packets
    data_sym_matrix = mapper_QPSK(data_bits_mat);
    
    % Create logical index for training symbols
    train_indices = mod(1:nDataOFDM, conf.train_period) == 0;
    
    % Total number of symbols including training symbols
    total_symbols = nDataOFDM + sum(train_indices);
    
    % Preallocate the all_ofdm_freq array
    all_ofdm_freq = zeros(N, total_symbols);
    
    % Fill in data symbols
    all_ofdm_freq(:, 1:nDataOFDM) = data_sym_matrix;
    
    % Insert training symbols at appropriate positions
    train_positions = cumsum(train_indices);
    all_ofdm_freq(:, train_positions + (1:nDataOFDM)) = trainSym(:, train_indices);
    conf.txSymbolsOFDM = all_ofdm_freq;
    
    % IFFT + CP
    s = ifft(all_ofdm_freq, [], 1); % Perform IFFT along the first dimension
    cp = s(end-Ncp+1:end, :);      % Extract cyclic prefix for all symbols
    txBlocks = [cp; s];     

    ofdm_bb = txBlocks(:);
    
    % Resampling to Audio
    ofdm_bb_os = ofdm_tx_resampling(ofdm_bb, conf);
    
    % SC Preamble
    preBits = preamble_generate(conf.sc.nsyms);
    preSyms = 2 * preBits - 1; % BPSK
    pre_up  = upsample(preSyms, conf.sc.os_factor);
    preamble_bb = conv(pre_up, conf.sc.txpulse);
    
    conf.preamble_bb = preamble_bb;
    
    % Normalization 
    % Average powers
    P_preamble = mean(abs(preamble_bb).^2);
    P_ofdm     = mean(abs(ofdm_bb_os).^2);
    
    % Scaling factors to have Average Power = 1
    scale_pre  = 1 / sqrt(P_preamble);
    scale_ofdm = 1 / sqrt(P_ofdm);
    
    % Apply normalization
    preamble_norm = preamble_bb * scale_pre;
    ofdm_norm     = ofdm_bb_os * scale_ofdm;
    
    % Save the normalized preamble in the config
    conf.sc.preamble_bb = preamble_norm;
    
    tx_baseband = [preamble_norm; ofdm_norm];
    
    % Upconversion
    n = (0:length(tx_baseband)-1).';
    audio_carrier = exp(1j * 2 * pi * fc * n / fs_audio);
    passband_signal = real(tx_baseband .* audio_carrier);
    
    % Final Scaling to avoid audio clipping
    % Normalize so that the maximum peak is 0.9
    max_val = max(abs(passband_signal));
    
    % Protection to avoid division by zero if the signal is null
    if max_val == 0
        txsignal = passband_signal;
    else
        txsignal = (passband_signal / max_val) * 0.9;
    end
end
