function [txsignal, conf] = tx_comb(txbits, conf)
    % TX with:
    % - 1 full OFDM training symbol
    % - OFDM data symbols with comb pilots in frequency

    % Basic constants
    M           = 2^conf.modulation_order;    % 4 for QPSK
    bitsPerSym  = conf.modulation_order;      % 2
    N           = conf.ofdm.ncarrier;
    Ncp         = conf.ofdm.cplen;
    fs_audio    = conf.f_s;
    fc          = conf.f_c;
    os_sc       = conf.sc.os_factor;
    
    %% Comb pilot pattern
    pilotSpacing       = 10;                  % pilot every 10 subcarriers
    pilotIdx           = 1:pilotSpacing:N;
    conf.ofdm.pilotIdx = pilotIdx;
    Npil               = numel(pilotIdx);
    
    dataIdx            = setdiff(1:N, pilotIdx);
    conf.ofdm.dataIdx  = dataIdx;
    Ndata              = numel(dataIdx);
    
    % pilot values (BPSK, all +1)
    pilotValues            = ones(Npil, 1);
    conf.ofdm.pilotValues  = pilotValues;
    
    %% Map data bits to QPSK
    qpskSymsData = mapper_QPSK(txbits);
    
    if mod(length(qpskSymsData), Ndata) ~= 0
        error('txofdm: number of QPSK symbols (%d) is not multiple of Ndata=%d', ...
              length(qpskSymsData), Ndata);
    end
    nDataOfdmSym = length(qpskSymsData) / Ndata;
    
    dataSymsMatrix = reshape(qpskSymsData, Ndata, nDataOfdmSym);
    
    %% Build OFDM data symbols with comb pilots
    ofdmDataSymbols = zeros(N, nDataOfdmSym);
    for s = 1:nDataOfdmSym
        symFreq            = zeros(N,1);
        symFreq(pilotIdx)  = pilotValues;
        symFreq(dataIdx)   = dataSymsMatrix(:,s);
        ofdmDataSymbols(:,s) = symFreq;
    end
    
    %% Full training OFDM symbol
    trainBits    = preamble_generate(N * bitsPerSym);
    trainSymFreq = mapper_QPSK(trainBits);
    
    conf.ofdm.trainSymFreq = trainSymFreq;
    
    % [training | data1 | data2 | ...]
    ofdmSymbols = [trainSymFreq, ofdmDataSymbols];
    nOfdmSymTot = size(ofdmSymbols, 2);
    conf.ofdm.txSym = ofdmSymbols;
    
    %% IFFT + CP insertion
    txBlocks = zeros(N+Ncp, nOfdmSymTot);
    for i = 1:nOfdmSymTot
        sTD = ifft(ofdmSymbols(:,i));
        cp  = sTD(end-Ncp+1:end);
        txBlocks(:,i) = [cp; sTD];
    end
    
    % serialize
    ofdm_bb = txBlocks(:);
    
    %% Resample to audio Fs
    ofdm_bb_os = ofdm_tx_resampling(ofdm_bb, conf);
    
    %% Single-carrier preamble
    preambleBits = preamble_generate(conf.sc.nsyms);
    bpskSyms     = 2*preambleBits - 1;
    
    bpsk_up     = upsample(bpskSyms, os_sc);
    preamble_bb = conv(bpsk_up, conf.sc.txpulse);
    
    %% Power normalization
    P_preamble = mean(abs(preamble_bb).^2);
    P_ofdm     = mean(abs(ofdm_bb_os).^2);
    
    scale_pre  = 1 / sqrt(P_preamble);
    scale_ofdm = 1 / sqrt(P_ofdm);
    
    preamble_norm = preamble_bb * scale_pre;
    ofdm_norm     = ofdm_bb_os * scale_ofdm;
    
    conf.sc.preamble_bb = preamble_norm;
    
    tx_baseband = [preamble_norm; ofdm_norm];
    
    %% Upconversion to passband
    n              = (0:length(tx_baseband)-1).';
    audio_carrier  = exp(1j*2*pi*fc*n/fs_audio);
    passband_signal = real(tx_baseband .* audio_carrier);
    
    % peak normalization
    max_val = max(abs(passband_signal));
    if max_val == 0
        txsignal = passband_signal;
    else
        txsignal = (passband_signal / max_val) * 0.9;
    end
end
