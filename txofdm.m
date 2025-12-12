function [txsignal, conf] = txofdm(txbits, conf)
    % TXOFDM - Transmitter function for OFDM modulation
    %
    %   [txsignal, conf] = txofdm(txbits, conf) implements a transmitter
    %   that can operate in different modes based on the configuration.
    %
    %   txbits : Input bits to be transmitted
    %   conf   : Configuration structure containing parameters for transmission
    %
    % Output:
    %   txsignal : The transmitted signal in the time domain
    %   conf     : Updated configuration structure with any modifications
    %
    % Modes:
    %   - task1: Standard OFDM transmission with a training symbol.
    %   - extratask1: Transmission of image / packets with periodic training OFDM.
    
    % Default transmitter mode if not specified
    if ~isfield(conf, 'tx_mode')
        conf.tx_mode = "task1"; % Set default mode to task1
    end

    switch string(conf.tx_mode)
        case "task1"
            [txsignal, conf] = tx_task1(txbits, conf);

        case "extratask1"
            [txsignal, conf] = tx_extratask1(txbits, conf);

        case "combtype"
            [txsignal, conf] = tx_comb_alt(txbits, conf); 

        otherwise
            error("txofdm: Unknown transmitter mode '%s'", conf.tx_mode); % Error for unknown mode
    end
end


%% ==================================================================
%  Trasmitter for most of the tasks
%  ==================================================================
function [txsignal, conf] = tx_task1(txbits,conf)
    % Digital Transmitter
    %
    %   [txsignal conf] = tx_task1(txbits,conf) implements a complete transmitter
    %   using a single carrier preamble followed by OFDM data in digital domain.
    %
    %   txbits  : Information bits
    %   conf    : Universal configuration structure
    
    % Useful constants
    M           = 2^conf.modulation_order;    
    bitsPerSym  = conf.modulation_order;      
    N           = conf.ofdm.ncarrier;         
    Ncp         = conf.ofdm.cplen;            
    fs_audio    = conf.f_s;                   
    fc          = conf.f_c;                   
    os_sc       = conf.sc.os_factor;
    
    % mapping
    qpskSyms = mapper_QPSK(txbits);
    
    % Symbols in OFDM matrix -> N x Nsym
    nOfdmSym = length(qpskSyms)/N;
    ofdmSymbols = reshape(qpskSyms,N,nOfdmSym); %512 x 50

    % Insert known training OFDM symbol
    trainBits = preamble_generate(N * bitsPerSym);
    trainSymFreq = mapper_QPSK(trainBits); 
    
    % Save in conf train symbol to use it in RX
    conf.ofdm.trainSymFreq = trainSymFreq;

    % Append training symbol to the OFDM symbols
    ofdmSymbols = [trainSymFreq, ofdmSymbols]; % Prepend training symbol
    conf.ofdm.txSym = ofdmSymbols;             % Used for task 3
    nOfdmSym = nOfdmSym + 1; % Update the number of OFDM symbols
    
    % IFFT + CP insertion
    s = ifft(ofdmSymbols, [], 1); % Perform IFFT along the first dimension
    cp = s(end-Ncp+1:end, :);      % Extract cyclic prefix for all symbols
    txBlocks = [cp; s];            % Combine CP and IFFT output
    
    % Serialize OFDM blocks
    ofdm_bb = txBlocks(:);
    
    % Resampling at audio frequency
    ofdm_bb_os = ofdm_tx_resampling(ofdm_bb, conf);     % now at 48 kHz
    %plot(ofdm_bb_os)
    
    % Generate single carrier BPSK preamble
    preambleBits = preamble_generate(conf.sc.nsyms);    % LFSR preamble

    % BPSK mapping: 0 -> +1, 1 -> -1
    bpskSyms = 2*preambleBits - 1;

    % Upsampling at audio frequency
    bpsk_up = upsample(bpskSyms, os_sc);

    % Pulse shaping with RRC pulse
    preamble_bb = conv(bpsk_up, conf.sc.txpulse);
    conf.sc.preamble_bb = preamble_bb;

    % NORMALIZATION 
    % Average powers
    P_preamble = mean(abs(preamble_bb).^2);
    P_ofdm     = mean(abs(ofdm_bb_os).^2);
    
    % Scaling factors to have Average Power = 1
    scale_pre  = 1 / sqrt(P_preamble);
    scale_ofdm = 1 / sqrt(P_ofdm);
    
    % Apply normalization
    preamble_norm = preamble_bb * scale_pre;
    ofdm_norm     = ofdm_bb_os * scale_ofdm;

    % IMPORTANT: Save the normalized preamble in the config
    conf.sc.preamble_bb = preamble_norm;

    tx_baseband = [preamble_norm; ofdm_norm];
    
    % Upconversion
    n = (0:length(tx_baseband)-1).';
    audio_carrier = exp(1j*2*pi*fc*n/fs_audio);
    passband_signal = real(tx_baseband .* audio_carrier);

    % Final Scaling to avoid audio clipping
    % Normalize so that the maximum peak is 0.9
    max_val = max(abs(passband_signal));
    
    % Protection to avoid division by zero if the signal is null
    if max_val == 0
        txsignal = passband_signal;
    else
        txsignal = (passband_signal / max_val) * 0.8;
    end
    
    % % Plot the transmitted signal
    % figure;
    % plot(txsignal);
    % title('Transmitted Signal');
    % xlabel('Sample Index');
    % ylabel('Amplitude');
    % grid on;
end

%% ==================================================================
%  Trasmitter for ExtraTask1: image / packets transmission 
%  ==================================================================
function [txsignal, conf] = tx_extratask1(bits_in, conf)
    % TX_EXTRATASK1 - Transmitter for image packets with periodic training OFDM
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
    train_period = conf.ofdm.train_period;

    % Training periodicity
    if ~isfield(conf.ofdm, "train_period")
        conf.ofdm.train_period = 20;   % every 20 data-OFDM symbols
    end
    
    % --- SCRAMBLER START ---
    % We use a specific seed (42) so the Receiver can reproduce the sequence.
    rng(42); 
    
    % Generate a random sequence of 0s and 1s the same length as your data
    scramble_seq = randi([0 1], length(bits_in), 1);
    conf.ofdm.scramble_seq = scramble_seq;
    
    % XOR your image bits with this random sequence
    bits_scrambled = xor(bits_in, scramble_seq);
    
    % Replace the original bits with the scrambled ones
    bits_in = double(bits_scrambled);

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
    
    [N, nData] = size(data_sym_matrix);

    % Determine where training symbols go
    train_flags = (mod(1:nData, train_period) == 1);    % logical positions
    nTrain = sum(train_flags);
    
    total_symbols = nData + nTrain;
    
    % ----- 2. Compute output positions for data -----
    % Shift each data column by number of training inserted before it
    data_shifts = cumsum(train_flags);                  % grows when training occurs
    data_pos = (1:nData) + data_shifts;                 % final output indices
    
    % ----- 3. Compute output positions for training symbols -----
    train_pos = data_pos(train_flags) - 1;              % training goes BEFORE each flagged data symbol
    
    % ----- 4. Build OFDM freq matrix -----
    all_ofdm_freq = zeros(N, total_symbols);
    
    % Fill data in vectorized manner
    all_ofdm_freq(:, data_pos) = data_sym_matrix;
    
    % Replicate training symbol to match number of trainings
    train_block = repmat(trainSym, 1, nTrain);
    
    % Fill training
    all_ofdm_freq(:, train_pos) = train_block;
    
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
    plot(txsignal)
end


%% ==================================================================
%  Transmitter for comb type training
% ===================================================================

function [txsignal, conf] = tx_comb_alt(txbits, conf)
    % TX with:
    % - 1 full OFDM training symbol (all subcarriers known)
    % - OFDM data symbols with diagonal comb pilots

    % Basic constants
    bitsPerSym = conf.modulation_order;   % 2 for QPSK
    N          = conf.ofdm.ncarrier;
    Ncp        = conf.ofdm.cplen;
    fs_audio   = conf.f_s;
    fc         = conf.f_c;
    os_sc      = conf.sc.os_factor;
    
    % ---------- Pilot pattern (diagonal comb) ----------
    pilotSpacing = conf.ofdm.pilotSpacing;
    conf.ofdm.pilotSpacing = pilotSpacing;
    
    baseIdx = 1:pilotSpacing:N;
    Npil    = numel(baseIdx);
    Ndata   = N - Npil;
    
    conf.ofdm.pilotBaseIdx = baseIdx;
    
    % shift per OFDM symbol (1 subcarrier = nice diagonal)
    pilotShift = 1;
    conf.ofdm.pilotShift = pilotShift;
    
    % pilot values (BPSK, all +1)
    pilotValues = ones(Npil, 1);
    conf.ofdm.pilotValues = pilotValues;
    
    % ---------- Map data bits to QPSK (WITH PADDING FOR Ndata) ----------
    bitsPerSym = conf.modulation_order;   % 2 for QPSK
    bitsPerOfdmDataSym = Ndata * bitsPerSym;
    
    nbits = numel(txbits);
    
    % pad bits so we have an integer number of OFDM data symbols
    nDataOfdmSym = ceil(nbits / bitsPerOfdmDataSym);
    nbits_pad = nDataOfdmSym * bitsPerOfdmDataSym - nbits;
    
    padBits = randi([0 1], nbits_pad, 1);
    txbits_pad = [txbits; padBits];
    conf.nbits_padded = numel(txbits_pad);
    conf.npad_bits_comb = nbits_pad;
    
    % map + reshape
    qpskSymsData = mapper_QPSK(txbits_pad);
    dataSymsMatrix = reshape(qpskSymsData, Ndata, nDataOfdmSym);

    % ---------- Build OFDM data symbols with diagonal pilots ----------
    ofdmDataSymbols = zeros(N, nDataOfdmSym);
    for s = 1:nDataOfdmSym
        shift      = (s-1)*pilotShift;
        pilotIdx_s = mod(baseIdx - 1 + shift, N) + 1;
        pilotIdx_s = sort(pilotIdx_s);
        
        dataIdx_s = setdiff(1:N, pilotIdx_s);
        if numel(dataIdx_s) ~= Ndata
            error('tx_comb: dataIdx_s size (%d) != Ndata (%d) in symbol %d.', ...
                  numel(dataIdx_s), Ndata, s);
        end
        
        symFreq = zeros(N, 1);
        symFreq(pilotIdx_s) = pilotValues;
        symFreq(dataIdx_s)  = dataSymsMatrix(:, s);
        
        ofdmDataSymbols(:, s) = symFreq;
    end
    
    % ---------- Full training OFDM symbol ----------
    trainBits    = preamble_generate(N * bitsPerSym);
    trainSymFreq = mapper_QPSK(trainBits);
    
    conf.ofdm.trainSymFreq = trainSymFreq;
    
    % [training | data1 | data2 | ...]
    ofdmSymbols = [trainSymFreq, ofdmDataSymbols];
    nOfdmSymTot = size(ofdmSymbols, 2);
    conf.ofdm.txSym = ofdmSymbols;
    
    % ---------- IFFT + CP ----------
    txBlocks = zeros(N+Ncp, nOfdmSymTot);
    for i = 1:nOfdmSymTot
        sTD = ifft(ofdmSymbols(:, i));
        cp  = sTD(end-Ncp+1:end);
        txBlocks(:, i) = [cp; sTD];
    end
    
    % Serialize OFDM blocks
    ofdm_bb = txBlocks(:);
    
    % ---------- Resample to audio Fs ----------
    ofdm_bb_os = ofdm_tx_resampling(ofdm_bb, conf);
    
    % ---------- Single-carrier preamble ----------
    preambleBits = preamble_generate(conf.sc.nsyms);
    bpskSyms     = 2 * preambleBits - 1;
    
    bpsk_up     = upsample(bpskSyms, os_sc);
    preamble_bb = conv(bpsk_up, conf.sc.txpulse);
    
    % ---------- Power normalization ----------
    P_preamble = mean(abs(preamble_bb).^2);
    P_ofdm     = mean(abs(ofdm_bb_os).^2);
    
    scale_pre  = 1 / sqrt(P_preamble);
    scale_ofdm = 1 / sqrt(P_ofdm);
    
    preamble_norm = preamble_bb * scale_pre;
    ofdm_norm     = ofdm_bb_os * scale_ofdm;
    
    conf.sc.preamble_bb = preamble_norm;
    
    tx_baseband = [preamble_norm; ofdm_norm];
    
    % ---------- Upconversion ----------
    n             = (0:length(tx_baseband)-1).';
    audio_carrier = exp(1j * 2*pi*fc*n/fs_audio);
    passband_signal = real(tx_baseband .* audio_carrier);
    
    max_val = max(abs(passband_signal));
    if max_val == 0
        txsignal = passband_signal;
    else
        txsignal = (passband_signal / max_val) * 0.9;
    end
end

