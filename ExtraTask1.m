close all; clear all; clc;

%% ===================== SYSTEM CONFIG ============================
conf.audiosystem = 'emulator';

conf.rx_mode = "task2";

%Emulator configuration
conf.emulator_idx = 2;
conf.emulator_snr = 20;

% General Parameters
conf.f_c   = 8000;

% Preamble
conf.sc.f_sym = 1000;
conf.sc.nsyms = 500;

% OFDM
conf.ofdm.bandwidth = 2000;
conf.ofdm.ncarrier  = 512;
conf.ofdm.cplen     = 256;
conf.modulation_order = 2;
conf.ofdm.OFDMnSyms = 50;

% Audio settings
conf.f_s = 48000;
conf.bitsps = 16;

% Calculations to be done once
conf.ofdm.spacing   = conf.ofdm.bandwidth/conf.ofdm.ncarrier;
conf.sc.os_factor   = conf.f_s/conf.sc.f_sym;

if mod(conf.sc.os_factor,1) ~= 0
   disp('WARNING: Sampling rate must be a multiple of the symbol rate for single carrier system'); 
end

conf.ofdm.os_factor = conf.f_s/(conf.ofdm.ncarrier*conf.ofdm.spacing);

% Pregenerate useful data
conf.sc.txpulse_length = 20*conf.sc.os_factor;
conf.sc.txpulse        = rrc(conf.sc.os_factor, 0.22, conf.sc.txpulse_length);

rng(0);


%% ===================== LOAD IMAGE ============================
img_tx = imread("img.png");   % grayscale image (uint8)
[H, W] = size(img_tx);

tx_bits = image_encoder(img_tx);
nbits = length(tx_bits);
conf.nbits = nbits;

% Add padding
bitsPerOFDM = conf.ofdm.ncarrier * conf.modulation_order;   % = 1024
needed = ceil(nbits / bitsPerOFDM) * bitsPerOFDM;

padding = needed - nbits;
tx_bits = [tx_bits ; zeros(padding,1)];

conf.nbits = length(tx_bits);

fprintf("Padding applied: %d bits → final %d bits\n", padding, conf.nbits);

fprintf("Image: %d x %d  →  %d bits\n", H, W, nbits);


%% ===================== TRANSMISSION ============================
disp('Start OFDM Transmission')

[rx_bits, conf] = audiotrans_packets(tx_bits, conf);

%% ===================== BER ================================
rx_bits = rx_bits(1:nbits);
tx_bits = tx_bits(:);

L = min(length(tx_bits), length(rx_bits));
BER = sum(tx_bits(1:L) ~= rx_bits(1:L)) / L;

fprintf("BER = %.6f\n", BER);

%% ===================== IMAGE DECODING =====================
img_rx = image_decoder(rx_bits(1:8*H*W), [H W]);

figure;
subplot(1,2,1); imshow(img_tx); title("TX Image");
subplot(1,2,2); imshow(img_rx); title("RX Image");
