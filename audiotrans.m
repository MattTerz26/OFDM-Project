% % % % %
% Wireless Receivers: algorithms and architectures
% Audio Transmission Framework 

clc;
clear all;
close all;

% Configuration Values

% Options for transmission are : 
% emulator: use a channel emulator with 5 different configurations.
% audio: use the loudspeaker and microphone for the data transmission
conf.audiosystem = 'emulator'; 

% Task reciever selection
conf.rx_mode = "task2";

%Emulator configuration
conf.emulator_idx = 2; % 1 to 5 yields different channels
conf.emulator_snr = 30;

% General parameters 
OFDMnSyms = 20;
conf.nbits   = 512*2*OFDMnSyms;  % number of bits 
conf.f_c     = 6000;

% Preamble
conf.sc.f_sym = 100;  % symbol rate
conf.sc.nsyms = 500;   % number of preamble symbols

%OFDM
conf.ofdm.bandwidth = 2000; %f_s/bw should be an integer for resampling
conf.ofdm.ncarrier  = 512;
conf.ofdm.cplen     = 256;
conf.modulation_order = 2; % 2 for QPSK

% Fix audio settings 
conf.f_s = 48000; % audio sampling rate, fixed by audiocard
conf.bitsps = 16;   % bits per audio sample

% all calculations that you only have to do once
conf.ofdm.spacing  = conf.ofdm.bandwidth/conf.ofdm.ncarrier;
conf.sc.os_factor  = conf.f_s/conf.sc.f_sym;

if mod(conf.sc.os_factor,1) ~= 0
   disp('WARNING: Sampling rate must be a multiple of the symbol rate for single carrier system'); 
end

conf.ofdm.os_factor = conf.f_s/(conf.ofdm.ncarrier*conf.ofdm.spacing);

% Pregenerate useful data
conf.sc.txpulse_length = 20*conf.sc.os_factor;
conf.sc.txpulse    = rrc(conf.sc.os_factor,0.22,conf.sc.txpulse_length);

disp('Start OFDM Transmission')

% Generate random data
rng(0);
txbits = randi([0 1],conf.nbits,1);

% Transmit Function
[txsignal, conf] = txofdm(txbits,conf);

rawtxsignal = [ zeros(conf.f_s,1) ; txsignal ; zeros(conf.f_s,1) ];
rawtxsignal = [  rawtxsignal  zeros(size(rawtxsignal)) ];

% MATLAB audio mode
switch(conf.audiosystem)
    case 'emulator'
        rxsignal = channel_emulator(rawtxsignal(:,1),conf);
    case 'audio'
        % % % % % % % % % % % % %
        % Begin
        % Audio Transmission    
       
        txdur       = length(rawtxsignal)/conf.f_s; % calculate length of transmitted signal
        audiowrite('out.wav',rawtxsignal,conf.f_s)
        
        disp('MATLAB generic');
        playobj = audioplayer(rawtxsignal,conf.f_s,conf.bitsps);
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
     
% Receive Function
[rxbits, conf]       = rxofdm(rxsignal,conf);

res.biterrors    = sum(rxbits ~= txbits);
ber = res.biterrors/length(rxbits)