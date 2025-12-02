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
    %trainBits = preamble_generate(N * bitsPerSym);
    %trainSymFreq = mapper_QPSK(trainBits); 
    trainSymFreq = (1+1j)/sqrt(2) * ones(N, 1);
    
    % Save in conf train symbol to use it in RX
    conf.ofdm.trainSymFreq = trainSymFreq;
    
    % Append training symbol to the OFDM symbols
    ofdmSymbols = [trainSymFreq, ofdmSymbols]; % Prepend training symbol
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

    %Mean Power 
    P_preamble = mean(abs(preamble_bb).^2);
    P_ofdm = mean(abs(ofdm_bb_os).^2);
    
    %Normalization
    P0 = 1; 
    alpha_pre  = sqrt(P0 / P_preamble);
    alpha_ofdm = sqrt(P0 / P_ofdm);
    preamble_bb_norm = alpha_pre  * preamble_bb;
    %conf.sc.preamble_bb = preamble_bb_norm;
    ofdm_bb_os_norm  = alpha_ofdm * ofdm_bb_os;

    % Concatenate preamble and OFDM basebande
    tx_baseband = [preamble_bb; ofdm_bb_os];
    
    % Mix the signal with the carrier frequency
    n = (0:length(tx_baseband)-1).';
    audio_carrier = exp(1j*2*pi*fc*n/fs_audio);
    passband = real(tx_baseband .* audio_carrier);
    txsignal = passband / (max(abs(passband)) + 1e-12) * 0.9;
end