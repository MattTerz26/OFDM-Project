clear all;
close all;
clc;

save_fig = 1;

%% Base configuration
conf.audiosystem  = 'emulator';
conf.f_c          = 4000;

conf.ofdm.bandwidth   = 2000;
conf.ofdm.ncarrier    = 512;     % N
conf.ofdm.cplen       = 256;
conf.modulation_order = 2;       % QPSK

conf.sc.f_sym  = 1000;
conf.sc.nsyms  = 500;

conf.f_s    = 48000;
conf.bitsps = 16;

%% Derived parameters
conf.ofdm.spacing   = conf.ofdm.bandwidth / conf.ofdm.ncarrier;
conf.sc.os_factor   = conf.f_s / conf.sc.f_sym;
conf.ofdm.os_factor = conf.f_s / (conf.ofdm.ncarrier * conf.ofdm.spacing);

conf.sc.txpulse_length = 20 * conf.sc.os_factor;
conf.sc.txpulse        = rrc(conf.sc.os_factor, 0.22, conf.sc.txpulse_length);

conf.rx_mode = "comb";   % diagonal pilot system

%% Fixed SNR
SNR_dB = 20;
conf.emulator_snr = SNR_dB;

%% Sweep parameters (pilot spacing)
emulator_ids     = 1:5;
pilotSpacing_vec = 20:10:(conf.ofdm.ncarrier * 0.75);

nChan = numel(emulator_ids);
nSpac = numel(pilotSpacing_vec);

BER_mat  = zeros(nChan, nSpac);
ERR_mat  = zeros(nChan, nSpac);
NbitsMat = zeros(nChan, nSpac);

%% Deterministic bit pool
rng(12345);
MAX_BITS = 1e6;
bitpool  = randi([0 1], MAX_BITS, 1);

%% Spacing + channel sweep
for sIdx = 1:nSpac
    
    pilotSpacing = pilotSpacing_vec(sIdx);
    conf.ofdm.pilotSpacing = pilotSpacing;

    N           = conf.ofdm.ncarrier;
    bitsPerSym  = conf.modulation_order;

    pilotIdx_A = 1:pilotSpacing:N;
    Npil       = numel(pilotIdx_A);
    Ndata      = N - Npil;

    nDataOfdmSym = 50;
    conf.nbits   = bitsPerSym * Ndata * nDataOfdmSym;

    if conf.nbits > MAX_BITS
        error('bitpool too small (MAX_BITS=%d, needed=%d).', MAX_BITS, conf.nbits);
    end

    fprintf('\n==== SPACING %d ====\n', pilotSpacing);

    for cIdx = 1:nChan
        if cIdx == 4
            continue;   % skip channel 4 
        end

        conf.emulator_idx = emulator_ids(cIdx);

        txbits = bitpool(1:conf.nbits);

        [txsignal, conf] = tx_comb_alt(txbits, conf);

        rawtxsignal = [zeros(conf.f_s,1); txsignal; zeros(conf.f_s,1)];
        rawtxsignal = [rawtxsignal zeros(size(rawtxsignal))];

        rxsignal = channel_emulator(rawtxsignal(:,1), conf);
        [rxbits, conf] = rx_comb_alt(rxsignal, conf);

        biterrors = sum(rxbits ~= txbits);
        BER       = biterrors / length(rxbits);

        BER_mat(cIdx, sIdx)  = BER;
        ERR_mat(cIdx, sIdx)  = biterrors;
        NbitsMat(cIdx, sIdx) = length(rxbits);

        fprintf('Ch %d -> BER %.3e  (%d/%d)\n', ...
            conf.emulator_idx, BER, biterrors, length(rxbits));
    end
end

%% Summary
fprintf('\n===== SUMMARY (SNR = %d dB) =====\n', SNR_dB);
for sIdx = 1:nSpac
    ps = pilotSpacing_vec(sIdx);
    fprintf('\nSpacing %d:\n', ps);
    for cIdx = 1:nChan
        if cIdx == 4
            continue;   % skip channel 4 
        end
        fprintf('Ch %d -> BER %.3e  (%d/%d)\n', ...
            emulator_ids(cIdx), BER_mat(cIdx, sIdx), ...
            ERR_mat(cIdx, sIdx), NbitsMat(cIdx, sIdx));
    end
end

%% Plot BER vs spacing
if save_fig
    figure; hold on;

    markers = {'o','s','d','^','v'};
    lines   = {'-','--','-.',':','-'};

    for cIdx = 1:nChan
        if cIdx == 4
            continue;
        end
        semilogy(pilotSpacing_vec, BER_mat(cIdx,:), ...
            lines{mod(cIdx-1,numel(lines))+1}, ...
            'Marker', markers{mod(cIdx-1,numel(markers))+1}, ...
            'LineWidth', 1.5, 'MarkerSize', 7);
    end

    grid on;
    xlabel('pilotSpacing');
    ylabel('BER');
    title(sprintf('BER vs pilotSpacing @ %d dB', SNR_dB));
    yline(1e-2, ':', 'BER=10^{-2}', 'HandleVisibility','off');

    legend(arrayfun(@(id) sprintf('Ch %d', id), emulator_ids, ...
        'UniformOutput', false), 'Location','southwest');

    saveas(gcf, 'ber_diag_pilot_spacing_lines.png');
end
