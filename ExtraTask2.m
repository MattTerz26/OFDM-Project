close all; clear all; clc;

%% 1. CONFIGURATION
rng(42);
conf.f_c     = 6000;
conf.f_s     = 48000;
conf.bitsps  = 16;
conf.Nsym    = 100;           
conf.nbits   = 512 * 2 * conf.Nsym;   
conf.enable_scrambler = true; % Toggle this ON

% SC Preamble
conf.sc.f_sym = 1000;          
conf.sc.nsyms = 500;     
conf.sc.os_factor = conf.f_s/conf.sc.f_sym;
conf.sc.txpulse_length = 20*conf.sc.os_factor;
conf.sc.txpulse = rrc(conf.sc.os_factor, 0.22, conf.sc.txpulse_length);

% OFDM
conf.ofdm.bandwidth = 2000;    
conf.ofdm.ncarrier  = 512;     
conf.ofdm.cplen     = 256;     
conf.modulation_order = 2;  

% all calculations that you only have to do once
conf.ofdm.spacing  = conf.ofdm.bandwidth/conf.ofdm.ncarrier;
conf.sc.os_factor  = conf.f_s/conf.sc.f_sym;
conf.ofdm.os_factor = conf.f_s/(conf.ofdm.ncarrier*conf.ofdm.spacing);

% Pregenerate useful data
conf.sc.txpulse_length = 20*conf.sc.os_factor;
conf.sc.txpulse    = rrc(conf.sc.os_factor,0.22,conf.sc.txpulse_length);

%% 2. GENERATE SIGNAL (TX)
txbits_in = preamble_generate(conf.nbits); 

% Generate Waveform
[txsignal, conf] = txofdm(txbits_in, conf); 

rawtxsignal = [ zeros(conf.f_s,1) ; txsignal ; zeros(conf.f_s,1) ];
rawtxsignal = [  rawtxsignal  zeros(size(rawtxsignal)) ];

%% 3. AUDIO TRANSMISSION 

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

%% 4. RX PROCESSING
try
    [rxbits, conf, H_est, Z_data] = rx_extratask2(rawrxsignal, conf);
    disp('RX Task 2 Completed Successfully.');
catch ME
    error('RX Failed: %s', ME.message);
end

%% 5. ANALYSIS & CHANNEL VISUALIZATION

% --- 1. BER Calculation ---
% Ensure we compare the same number of bits
nBitsComp = min(length(rxbits), length(txbits_in));
ber = sum(rxbits(1:nBitsComp) ~= txbits_in(1:nBitsComp)) / nBitsComp;
fprintf('\n<strong>BER = %.5f</strong>\n', ber);

% --- 2. Reconstruct Instantaneous Channel (Decision Directed) ---
% To see evolution, we need: H_inst = Z_received / X_transmitted
% We regenerate the X_transmitted from the original bits
bitsPerSym = conf.modulation_order;
N = conf.ofdm.ncarrier;

% Reshape original bits to symbols
if mod(length(txbits_in), bitsPerSym) ~= 0
    % Handle case where padding might be needed (unlikely if config is correct)
    txbits_in = [txbits_in; zeros(mod(length(txbits_in), bitsPerSym), 1)];
end

% Map bits to QPSK Symbols
tx_syms_serial = mapper_QPSK(txbits_in);

% Reshape to Grid [N x nDataSym]
% Ensure dimensions match Z_data (Received Data Grid)
nDataSym = size(Z_data, 2);
% Truncate or pad TX symbols to match RX length exactly
if length(tx_syms_serial)/N >= nDataSym
    tx_syms_grid = reshape(tx_syms_serial(1:N*nDataSym), N, nDataSym);
    
    % Compute Instantaneous Channel Evolution
    H_evolution = Z_data ./ tx_syms_grid;
else
    warning('Not enough TX bits to estimate full channel evolution.');
    H_evolution = Z_data; % Fallback
end

% --- 3. PLOTS ---

% A. Power Delay Profile (PDP) - To estimate Channel Length
% We use the clean estimate from the Preamble (H_est)
h_impulse = ifft(H_est, N, 1);       % Transform to time domain
pdp = abs(h_impulse).^2;
pdp_dB = 10*log10(pdp);
pdp_dB = pdp_dB - max(pdp_dB);       % Normalize to 0dB

figure;
plot(0:N-1, pdp_dB, 'LineWidth', 2);
yline(-20, 'r--', 'Threshold (-20dB)');
xlim([0 50]);                        % Zoom in on the first 50 taps (usually where audio reflections are)
title('Power Delay Profile (PDP)');
xlabel('Delay (Samples)');
ylabel('Normalized Power [dB]');
grid on;

% B. Channel Frequency Response (Snapshot)
H_dB = 20*log10(abs(fftshift(H_est)));
fAxis = linspace(-conf.ofdm.bandwidth/2, conf.ofdm.bandwidth/2, N);

figure;
plot(fAxis, H_dB, 'LineWidth', 1.5);
title('Channel Frequency Response (at Preamble)');
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
grid on;
xlim([min(fAxis) max(fAxis)]);

% C. Channel Evolution (Time-Frequency Heatmap)
% This shows how fading changes during the packet
figure;
imagesc(1:nDataSym, fAxis, 20*log10(abs(fftshift(H_evolution, 1))));
set(gca, 'YDir', 'normal'); % Correct Y-axis direction
colorbar;
clim([-20 5]); % Clamp colors to visualize fading clearly
title('Channel Evolution Magnitude [dB]');
xlabel('OFDM Symbol Index (Time)');
ylabel('Frequency [Hz]');

% D. Phase Evolution (Carrier Frequency Offset Check)
% We pick a central subcarrier and track its phase over time
subcarrier_idx = floor(N/2) + 10; % Pick a freq bin slightly off DC
phase_trace = unwrap(angle(H_evolution(subcarrier_idx, :)));

figure;
plot(1:nDataSym, rad2deg(phase_trace), '-o', 'MarkerSize', 3);
title(sprintf('Phase Evolution at Subcarrier %d', subcarrier_idx));
xlabel('OFDM Symbol Index');
ylabel('Phase [Degrees]');
grid on;
% Note: A constant slope here indicates Residual CFO.

% E. Constellation
figure;
plot(Z_data(:), '.', 'MarkerSize', 4);
hold on;
% Plot ideal QPSK points for reference
qpsk_ref = [1+1j, 1-1j, -1+1j, -1-1j] / sqrt(2); 
% Note: You might need to scale qpsk_ref based on your specific normalization
plot(qpsk_ref, 'rx', 'LineWidth', 2, 'MarkerSize', 10);
title(sprintf('Received Constellation (BER: %.5f)', ber));
axis equal; grid on;
hold off;