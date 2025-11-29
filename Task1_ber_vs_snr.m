clear all;
close all;
clc;

% ----------------- CONFIGURAZIONE BASE -----------------
conf.audiosystem    = 'emulator';
conf.emulator_idx   = 1;      % canale semplice per Task 1
conf.f_c            = 8000;

% General parameters 
conf.nbits   = 512*2*50;      % deve essere consistente con OFDM
conf.ofdm.bandwidth = 2000;
conf.ofdm.ncarrier  = 512;
conf.ofdm.cplen     = 256;
conf.modulation_order = 2;

% Preamble
conf.sc.f_sym = 1000;
conf.sc.nsyms = 500;

% Audio
conf.f_s    = 48000;
conf.bitsps = 16;

% Calcoli derivati
conf.ofdm.spacing   = conf.ofdm.bandwidth/conf.ofdm.ncarrier;
conf.sc.os_factor   = conf.f_s/conf.sc.f_sym;
conf.ofdm.os_factor = conf.f_s/(conf.ofdm.ncarrier*conf.ofdm.spacing);

conf.sc.txpulse_length = 20*conf.sc.os_factor;
conf.sc.txpulse        = rrc(conf.sc.os_factor,0.22,conf.sc.txpulse_length);

% Receiver mode (quello semplice per Task 1)
conf.rx_mode = "task1_basic";

% ----------------- VETTORE DI SNR -----------------
SNR_dB_vec = 0:1:40;      % ad esempio da 0 a 40 dB a passi di 5
BER_vec    = zeros(size(SNR_dB_vec));

% Per avere risultati riproducibili
rng(0);

for k = 1:length(SNR_dB_vec)
    conf.emulator_snr = SNR_dB_vec(k);

    % ----- Genera dati casuali -----
    txbits = randi([0 1], conf.nbits, 1);

    % ----- Trasmissione -----
    [txsignal, conf] = txofdm(txbits, conf);

    % Padding con zeri come nel template originale
    rawtxsignal = [ zeros(conf.f_s,1) ; txsignal ; zeros(conf.f_s,1) ];
    rawtxsignal = [ rawtxsignal  zeros(size(rawtxsignal)) ];

    % Canale emulato
    rxsignal = channel_emulator(rawtxsignal(:,1), conf);

    % ----- Ricezione -----
    [rxbits, conf] = rxofdm(rxsignal, conf);

    % ----- Calcolo BER -----
    biterrors = sum(rxbits ~= txbits);
    BER_vec(k) = biterrors/length(rxbits);

    fprintf('SNR = %2d dB  ->  BER = %.3e (errors = %d)\n', ...
             SNR_dB_vec(k), BER_vec(k), biterrors);
end

% ----------------- PLOT BER VS SNR -----------------
figure;
semilogy(SNR_dB_vec, BER_vec, '-o','LineWidth',1.5);
yline(0.01, '-w');
grid on;
xlabel('SNR [dB]');
ylabel('BER');
title('BER vs SNR (emulator idx 1)');
% save picture
saveas(gcf, 'ber_vs_snr_task1.png');