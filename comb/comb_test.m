clear all;
close all;
clc;

save_fig = 0;

%% Base configuration
conf.audiosystem  = 'emulator';
conf.emulator_idx = 5;
conf.f_c          = 8000;

%% OFDM parameters
conf.ofdm.bandwidth   = 2000;
conf.ofdm.ncarrier    = 512;   % N
conf.ofdm.cplen       = 256;
conf.modulation_order = 2;     % QPSK

%% Preamble
conf.sc.f_sym = 1000;
conf.sc.nsyms = 500;

%% Audio settings
conf.f_s    = 48000;
conf.bitsps = 16;

%% Derived parameters
conf.ofdm.spacing   = conf.ofdm.bandwidth / conf.ofdm.ncarrier;
conf.sc.os_factor   = conf.f_s / conf.sc.f_sym;
conf.ofdm.os_factor = conf.f_s / (conf.ofdm.ncarrier * conf.ofdm.spacing);

conf.sc.txpulse_length = 20 * conf.sc.os_factor;
conf.sc.txpulse        = rrc(conf.sc.os_factor, 0.22, conf.sc.txpulse_length);

%% Comb pilot setup
N          = conf.ofdm.ncarrier;
bitsPerSym = conf.modulation_order;

pilotSpacing = 10;
pilotIdx     = 1:pilotSpacing:N;
Npil         = numel(pilotIdx);
Ndata        = N - Npil;

nDataOfdmSym = 50;

conf.nbits = bitsPerSym * Ndata * nDataOfdmSym;

fprintf('N=%d, Npil=%d, Ndata=%d, nDataOfdmSym=%d\n', ...
        N, Npil, Ndata, nDataOfdmSym);
fprintf('conf.nbits = %d bits\n', conf.nbits);

conf.rx_mode = "comb";

%% Single SNR test
SNR_dB = 20;
conf.emulator_snr = SNR_dB;

rng(0);  % reproducible bits

%% Bit generation
txbits = randi([0 1], conf.nbits, 1);

%% Transmit
[txsignal, conf] = tx_comb(txbits, conf);

rawtxsignal = [zeros(conf.f_s,1); txsignal; zeros(conf.f_s,1)];
rawtxsignal = [rawtxsignal zeros(size(rawtxsignal))];

%% Channel
rxsignal = channel_emulator(rawtxsignal(:,1), conf);

%% Receive
[rxbits, conf] = rx_comb(rxsignal, conf);

%% BER
biterrors = sum(rxbits ~= txbits);
BER = biterrors / length(rxbits);

fprintf('SNR = %d dB -> BER = %.3e (errors = %d)\n', ...
        SNR_dB, BER, biterrors);

%% Optional plot
if save_fig
    figure;
    semilogy(SNR_dB, BER, 'o', 'LineWidth', 1.5);
    yline(1e-2, '-');
    grid on;
    xlabel('SNR [dB]');
    ylabel('BER');
    title(sprintf('BER at %d dB (Channel ID %d)', ...
          SNR_dB, conf.emulator_idx));
    saveas(gcf, 'ber_snr20_comb.png');
end
