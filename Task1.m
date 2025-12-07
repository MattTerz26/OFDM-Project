clear all;
close all;
clc;
save  = 0;

% ----------------- BASE CONFIGURATION -----------------
conf.audiosystem    = 'emulator';
conf.emulator_idx   = 1;      % simple channel for Task 1
conf.f_c            = 8000;

% General parameters 
conf.nbits   = 512*2*50;      % must be consistent with OFDM
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

% Derived calculations
conf.ofdm.spacing   = conf.ofdm.bandwidth/conf.ofdm.ncarrier;
conf.sc.os_factor   = conf.f_s/conf.sc.f_sym;
conf.ofdm.os_factor = conf.f_s/(conf.ofdm.ncarrier*conf.ofdm.spacing);

conf.sc.txpulse_length = 20*conf.sc.os_factor;
conf.sc.txpulse        = rrc(conf.sc.os_factor,0.22,conf.sc.txpulse_length);

% Receiver mode (the simple one for Task 1)
conf.rx_mode = "task1";

% ----------------- SNR VECTOR -----------------
SNR_dB_vec = 0:5:40;      % for example from 0 to 40 dB in steps of 5
BER_vec    = zeros(size(SNR_dB_vec));

% For reproducible results
rng(0);

for k = 1:length(SNR_dB_vec)
    conf.emulator_snr = SNR_dB_vec(k);

    % ----- Generate random data -----
    txbits = randi([0 1], conf.nbits, 1);

    % ----- Transmission -----
    [txsignal, conf] = txofdm(txbits, conf);

    % Padding with zeros as in the original template
    rawtxsignal = [ zeros(conf.f_s,1) ; txsignal ; zeros(conf.f_s,1) ];
    rawtxsignal = [ rawtxsignal  zeros(size(rawtxsignal)) ];

    % Emulated channel
    rxsignal = channel_emulator(rawtxsignal(:,1), conf);

    % ----- Reception -----
    [rxbits, conf] = rxofdm(rxsignal, conf);

    % ----- BER Calculation -----
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
title('BER vs SNR - Channel ID 1');

% save picture
if save
    saveas(gcf, 'ber_vs_snr_task1.png');
end
