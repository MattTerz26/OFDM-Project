function [rxbits, conf] = rxofdm(rxsignal, conf)
% Digital Receiver (dispatcher)
%
%   [rxbits conf] = rxofdm(rxsignal,conf) implements different receiver
%   variants depending on conf.rx_mode.
%
%   rxsignal : received real passband signal (from emulator/audio)
%   conf     : configuration structure
%   rxbits   : recovered information bits

    % Default receiver mode if not specified
    if ~isfield(conf, 'rx_mode')
        conf.rx_mode = "task1_basic";
    end

    switch string(conf.rx_mode)
        case "task1_basic"
            rxbits = rx_task1_basic(rxsignal, conf);

        otherwise
            error("rxofdm: Unknown receiver mode '%s'", conf.rx_mode);
    end
end


%% ========================================================================
%  Receiver for Task 1: simplest OFDM receiver (no phase tracking)
% =======================================================================
function rxbits = rx_task1_basic(rxsignal, conf)

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

    % Remove leading/trailing zeros
    rx_trim = rxsignal(fs_audio+1 : end-fs_audio);   % keep only txsignal part

    % Downconversion to complex baseband
    n = (0:length(rx_trim)-1).';
    carrier_rx = exp(-1j*2*pi*fc*n/fs_audio);       % minus sign for RX
    bb_rx = rx_trim(:) .* carrier_rx;

    % Low-pass filtering (extract baseband)
    bb_rx_filt = ofdmlowpass(bb_rx, conf, conf.ofdm.bandwidth);

    % Remove preamble in time-domain by knowing its length
    %    Preamble length:
    %       - conf.sc.nsyms symbols
    %       - upsampled by os_sc
    %       - convolved with RRC of length Lp = 2*txpulse_length + 1
    %       => Np = nsyms*os_sc + Lp - 1

    Lp = length(conf.sc.txpulse);
    preamble_len = conf.sc.nsyms*os_sc + Lp - 1;

    if length(bb_rx_filt) <= preamble_len
        error('rx_task1_basic: received signal too short compared to preamble length.');
    end

    ofdm_bb_os_rx = bb_rx_filt(preamble_len+1:end);

    % Resample from audio Fs to OFDM Fs
    ofdm_bb_rx = ofdm_rx_resampling(ofdm_bb_os_rx, conf);

    % Shorten to expected length (multiple of (N+Ncp)*nOfdmSymTot )
    total_ofdm_len = nOfdmSymTot * (N + Ncp);

    if length(ofdm_bb_rx) < total_ofdm_len
        error('rx_task1_basic: ofdm_bb_rx too short (got %d, need %d).', ...
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
        error('rx_task1_basic: not enough bits decoded (got %d, need %d).', ...
               length(rxbits_all), conf.nbits);
    end
    rxbits = rxbits(1:conf.nbits);

end
