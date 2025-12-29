% sweep_center_frequency_ber.m
% Sweep der Trägerfrequenz (TX task1, RX task2) von 2 kHz bis 10 kHz in
% 2-kHz-Schritten, BER-Auswertung und Plot-Speicherung.

clear; close all; clc;

%% Grundkonfiguration
conf.audiosystem    = 'audio';   % 'emulator' oder 'audio'
conf.emulator_idx   = 2;            % Kanal-ID für Emulator
conf.emulator_snr   = 20;           % SNR im Emulator [dB]
conf.tx_mode        = "task1";
conf.rx_mode        = "task2";
conf.ofdm.alpha     = 0.1;          % IIR-Faktor für Phase-Tracking (rx task2)

% Allgemeine Parameter
conf.nbits          = 512*2*40;     % muss Vielfaches von (ncarrier*modorder) sein
conf.ofdm.bandwidth = 2000;
conf.ofdm.ncarrier  = 512;
conf.ofdm.cplen     = 256;
conf.modulation_order = 2;          % QPSK

% Preamble
conf.sc.f_sym = 1000;
conf.sc.nsyms = 500;

% Audio
conf.f_s    = 48000;
conf.bitsps = 16;

% Abgeleitete Größen
conf.ofdm.spacing   = conf.ofdm.bandwidth/conf.ofdm.ncarrier;
conf.sc.os_factor   = conf.f_s/conf.sc.f_sym;
conf.ofdm.os_factor = conf.f_s/(conf.ofdm.ncarrier*conf.ofdm.spacing);

conf.sc.txpulse_length = 20*conf.sc.os_factor;
conf.sc.txpulse        = rrc(conf.sc.os_factor,0.22,conf.sc.txpulse_length);

% Frequenz-Sweep (kHz -> Hz)
fc_vec = 2000:2000:10000;  % [Hz]
nFc    = numel(fc_vec);
BER_vec = zeros(1, nFc);

% Gleiche Bits für alle Sweeps
rng(0);
txbits_ref = randi([0 1], conf.nbits, 1);

for k = 1:nFc
    conf_k = conf;
    conf_k.f_c = fc_vec(k);

    % TX
    [txsignal, conf_tx] = txofdm(txbits_ref, conf_k);

    % Padding (wie audiotrans)
    rawtxsignal = [ zeros(conf_tx.f_s,1) ; txsignal ; zeros(conf_tx.f_s,1) ];
    rawtxsignal = [ rawtxsignal  zeros(size(rawtxsignal)) ]; % 2-Kanal-Format

    % Kanal
    switch conf_tx.audiosystem
        case 'emulator'
            rxsignal = channel_emulator(rawtxsignal(:,1), conf_tx);
        case 'audio'
            playobj = audioplayer(rawtxsignal, conf_tx.f_s, conf_tx.bitsps);
            recobj  = audiorecorder(conf_tx.f_s, conf_tx.bitsps, 1);
            record(recobj);
            pause(2);
            playblocking(playobj);
            pause(2);
            stop(recobj);
            rawrxsignal  = getaudiodata(recobj,'int16');
            rxsignal = double(rawrxsignal(1:end))/double(intmax('int16'));
        otherwise
            error('Unknown audiosystem: %s', conf_tx.audiosystem);
    end

    % RX (task2)
    conf_rx = conf_tx;
    conf_rx.rx_mode = "task2";
    [rxbits, conf_rx] = rxofdm(rxsignal, conf_rx);

    % BER
    biterrors = sum(rxbits ~= txbits_ref);
    BER_vec(k) = biterrors/length(rxbits);
    fprintf('f_c = %5d Hz -> BER = %.3e (errors = %d)\n', fc_vec(k), BER_vec(k), biterrors);
end

% Bester Wert
[minBer, idxBest] = min(BER_vec);
bestFc = fc_vec(idxBest);
fprintf('\nBestes Ergebnis: f_c = %d Hz mit BER = %.3e\n', bestFc, minBer);

%% Plot BER vs Center Frequency
fig = figure('Name','BER vs Center Frequency');
semilogy(fc_vec, BER_vec, '-o','LineWidth',1.5);
grid on;
xlabel('Center frequency f_c [Hz]');
ylabel('BER');
title('BER vs f_c (TX task1, RX task2)');

% Plot speichern
outdir = fullfile(pwd, 'plots_frequency_sweep');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
saveas(fig, fullfile(outdir, sprintf('ber_vs_fc_%s.png', timestamp)));

fprintf('Plot gespeichert in %s\n', outdir);
