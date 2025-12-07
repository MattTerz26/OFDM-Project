function [txsignal, conf] = txofdm(txbits,conf)
    % Digital Transmitter
    %
    %   [txsignal conf] = tx(txbits,conf) implements a complete transmitter
    %   using a single carrier preamble followed by OFDM data in digital domain.
    %
    %   txbits  : Information bits
    %   conf    : Universal configuration structure
    
    % Useful constants
    M           = 2^conf.modulation_order;    % 4 for QPSK
    bitsPerSym  = conf.modulation_order;      % 2
    N           = conf.ofdm.ncarrier;         % 512
    Ncp         = conf.ofdm.cplen;            % 256
    fs_audio    = conf.f_s;                   % 48000 Hz
    fc          = conf.f_c;                   % 8000 Hz
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
    txBlocks = zeros(N+Ncp,nOfdmSym);
    for i = 1:nOfdmSym
        %IFFT
        s = ifft(ofdmSymbols(:,i));
        cp = s(end-Ncp+1:end);
        %combine payload and cp
        txBlocks(:,i) = [cp; s];
    end
    
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
    %plot(preamble_bb)

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
        txsignal = (passband_signal / max_val) * 0.9;
    end

end