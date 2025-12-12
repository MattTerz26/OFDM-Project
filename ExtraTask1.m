close all; clear all; clc;

%% ===================== SYSTEM CONFIG ============================
packet_version = 2;

% Emulator configuration
conf.audiosystem = 'audio';
conf.emulator_idx = 3;
conf.emulator_snr = 15;

% General Parameters
conf.f_c   = 6000;

% Preamble
conf.sc.f_sym = 1000;
conf.sc.nsyms = 500;

% OFDM
conf.ofdm.bandwidth = 2000;
conf.ofdm.ncarrier  = 512;
conf.ofdm.cplen     = 256;
conf.modulation_order = 2;
conf.ofdm.OFDMnSyms = 20;       % How many OFDM symbols per packet, packet version 1
conf.ofdm.train_period = 5;    % How many OFDM symbols per packet, packet version 2

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
img_tx = imread("img4.png");   % grayscale image (uint8)
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

switch packet_version
    case 1
        disp('Packet version 1 selected: ( p + training + data ) for N times');
        conf.rx_mode = "task2";         % Uses old reciever and transmitter 
        [rx_bits, conf] = audiotrans_packets(tx_bits, conf);
    case 2
        disp('Packet version 2 selected: (p + train + data1 + train + data2 ...)');
        conf.rx_mode = "extratask1";
        conf.tx_mode = "extratask1";
        [txsignal, conf] = txofdm(tx_bits, conf);
        rawtx = [zeros(conf.f_s,1); txsignal; zeros(conf.f_s,1)];
        rawtx = [rawtx zeros(size(rawtx))];

        % --- CHANNEL ---
        switch(conf.audiosystem)
            case 'emulator'
                rxsignal = channel_emulator(rawtx(:,1),conf);
            case 'audio'
                % % % % % % % % % % % % %
                % Begin
                % Audio Transmission    
               
                txdur       = length(rawtx)/conf.f_s; % calculate length of transmitted signal
                audiowrite('out.wav',rawtx,conf.f_s)
                
                disp('MATLAB generic');
                playobj = audioplayer(rawtx,conf.f_s,conf.bitsps);
                recobj  = audiorecorder(conf.f_s,conf.bitsps,1);
                record(recobj);
                pause(2);
                disp('Recording...');
                playblocking(playobj)
                pause(2);
                stop(recobj);
                disp('Recording ended')
                rawrxsignal  = getaudiodata(recobj,'int16');
                rawrxsignal     = double(rawrxsignal(1:end))/double(intmax('int16')) ;
                rxsignal = rawrxsignal; 
        
                %
                % End
                % Audio Transmission   
                % % % % % % % % % % % %
        end

        % --- RX ---
        [rx_bits, conf] = rxofdm(rxsignal, conf);         % rx_mode extratask1
    otherwise   
        disp('Warning: Unsupported packet version. Please check the configuration.');
end


%% ===================== BER ================================
rx_bits = rx_bits(1:nbits);
tx_bits = tx_bits(:);

L = min(length(tx_bits), length(rx_bits));
BER = sum(tx_bits(1:L) ~= rx_bits(1:L)) / L;

fprintf("BER = %.6f\n", BER);

%% ===================== IMAGE DECODING =====================
img_rx = image_decoder(rx_bits, [H W]);

figure;
subplot(1,2,1); imshow(img_tx); title("TX Image");
subplot(1,2,2); imshow(img_rx); title("RX Image");
