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
    
    % ---------- Map data bits to QPSK ----------
    qpskSymsData = mapper_QPSK(txbits);
    nDataOfdmSym = conf.nbits / (bitsPerSym * Ndata);
    
    if mod(length(qpskSymsData), Ndata) ~= 0
        error('tx_comb: number of QPSK symbols (%d) not multiple of Ndata=%d', ...
              length(qpskSymsData), Ndata);
    end
    
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
