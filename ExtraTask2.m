close all; clear all; clc;

%% GLOBAL PLOT SETTINGS (bigger & report-like)
set(groot, ...
    'defaultFigureColor', 'w', ...
    'defaultAxesFontSize', 18, ...
    'defaultAxesLineWidth', 1.6, ...
    'defaultLineLineWidth', 2.4, ...
    'defaultLineMarkerSize', 9, ...
    'defaultTextFontSize', 18, ...
    'defaultAxesTitleFontSizeMultiplier', 1.25, ...
    'defaultAxesLabelFontSizeMultiplier', 1.15);

% Figure scaling (bigger window by default)
set(groot, 'defaultFigureUnits', 'pixels');
set(groot, 'defaultFigurePosition', [120 80 1100 720]); % [x y width height]

%% 1. CONFIGURATION
rng(42);
conf.f_c     = 6000;
conf.f_s     = 48000;
conf.bitsps  = 16;
conf.Nsym    = 50;
conf.nbits   = 512 * 2 * conf.Nsym;
conf.enable_scrambler = true; % Toggle this ON

% SC Preamble
conf.sc.f_sym = 500;
conf.sc.nsyms = 100;
conf.sc.os_factor = conf.f_s/conf.sc.f_sym;
conf.sc.txpulse_length = 20*conf.sc.os_factor;
conf.sc.txpulse = rrc(conf.sc.os_factor, 0.22, conf.sc.txpulse_length);

% OFDM
conf.ofdm.bandwidth = 2000;
conf.ofdm.ncarrier  = 512;
conf.ofdm.cplen     = 256;
conf.modulation_order = 2;

% all calculations that you only have to do once
conf.ofdm.spacing   = conf.ofdm.bandwidth/conf.ofdm.ncarrier;
conf.sc.os_factor   = conf.f_s/conf.sc.f_sym;
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
rawrxsignal  = double(rawrxsignal(1:end))/double(intmax('int16')) ;
rxsignal = rawrxsignal;

%% 4. RX PROCESSING
try
    [rxbits, conf, H_est, Z_data, Y] = rx_extratask2(rawrxsignal, conf);
    disp('RX Task 2 Completed Successfully.');
catch ME
    error('RX Failed: %s', ME.message);
end

%% 5. ANALYSIS & CHANNEL VISUALIZATION

% --- 1. BER Calculation ---
nBitsComp = min(length(rxbits), length(txbits_in));
ber = sum(rxbits(1:nBitsComp) ~= txbits_in(1:nBitsComp)) / nBitsComp;
fprintf('\n<strong>BER = %.5f</strong>\n', ber);

% --- 2. Reconstruct Instantaneous Channel (Decision Directed) ---
bitsPerSym = conf.modulation_order;
N = conf.ofdm.ncarrier;

if mod(length(txbits_in), bitsPerSym) ~= 0
    txbits_in = [txbits_in; zeros(mod(length(txbits_in), bitsPerSym), 1)];
end

tx_syms_serial = mapper_QPSK(txbits_in);

nDataSym = size(Z_data, 2);
if length(tx_syms_serial)/N >= nDataSym
    tx_syms_grid = reshape(tx_syms_serial(1:N*nDataSym), N, nDataSym);
    H_evolution = Z_data ./ tx_syms_grid;
else
    warning('Not enough TX bits to estimate full channel evolution.');
    H_evolution = Z_data; % Fallback
end

% --- 3. PLOTS ---

% A. Power Delay Profile (PDP)
h_impulse = ifft(H_est, N, 1);
pdp = abs(h_impulse).^2;
pdp_dB = 10*log10(pdp);
pdp_dB = pdp_dB - max(pdp_dB);

figure;
plot(0:N-1, pdp_dB, 'LineWidth', 2);
yline(-20, 'r--', 'Threshold (-20dB)');
xlim([0 50]);
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
figure;
imagesc(1:nDataSym, fAxis, 20*log10(abs(fftshift(H_evolution, 1))));
set(gca, 'YDir', 'normal');
colorbar;
clim([-20 5]);
title('Channel Evolution Magnitude [dB]');
xlabel('OFDM Symbol Index (Time)');
ylabel('Frequency [Hz]');

% D. Phase Evolution
subcarrier_idx = floor(N/2) + 10;
phase_trace = unwrap(angle(H_evolution(subcarrier_idx, :)));

figure;
plot(1:nDataSym, rad2deg(phase_trace), '-o', 'MarkerSize', 3);
title(sprintf('Phase Evolution at Subcarrier %d', subcarrier_idx));
xlabel('OFDM Symbol Index');
ylabel('Phase [Degrees]');
grid on;

% E. Constellation
figure;
plot(Y(:), '.', 'MarkerSize', 6);
hold on;
qpsk_ref = [1+1j, 1-1j, -1+1j, -1-1j] / sqrt(2);
plot(qpsk_ref, 'rx', 'LineWidth', 2, 'MarkerSize', 12);
title(sprintf('Received Constellation (BER: %.5f)', ber));
axis equal; grid on;
hold off;

%% 6. SAVE FIGURES (dated subfolder)
baseDir = 'plot_extratask2';
if ~exist(baseDir, 'dir')
    mkdir(baseDir);
end

runStamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');  % e.g. 2025-12-15_22-13-05
saveDir  = fullfile(baseDir, runStamp);
mkdir(saveDir);

figHandles = findall(groot, 'Type', 'figure');

% Save with consistent size & good resolution
for i = 1:length(figHandles)
    f = figHandles(i);

    % Ensure figure uses the larger default size when saving
    set(f, 'Units', 'pixels');
    pos = get(f, 'Position');
    if pos(3) < 1000 || pos(4) < 650
        set(f, 'Position', [pos(1) pos(2) 1100 720]);
    end

    outName = fullfile(saveDir, sprintf('figure_%02d.png', i));
    exportgraphics(f, outName, 'Resolution', 200);  % nicer than saveas
end

fprintf('Saved %d figures to: %s\n', length(figHandles), saveDir);
