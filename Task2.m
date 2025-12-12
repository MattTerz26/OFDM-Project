clear all;
close all;
clc;

save = 0; % Set this variable to control saving of plots

%% ----------------- BASE CONFIGURATION -----------------
conf.audiosystem    = 'emulator';
conf.emulator_idx   = 2;        % Channel ID 2 
conf.emulator_snr   = 10;      % Let's start with 100 dB

% General parameters 
conf.nbits   = 500*2*50;       
conf.f_c     = 8000;           

% Single-carrier preamble
conf.sc.f_sym = 1000;          
conf.sc.nsyms = 500;           

% OFDM
conf.ofdm.bandwidth = 2000;    
conf.ofdm.ncarrier  = 500;     
conf.ofdm.cplen     = 250;     
conf.modulation_order = 2;     

% Audio
conf.f_s    = 48000;           
conf.bitsps = 16;              

% Derived
conf.ofdm.spacing   = conf.ofdm.bandwidth/conf.ofdm.ncarrier;
conf.sc.os_factor   = conf.f_s/conf.sc.f_sym;
conf.ofdm.os_factor = conf.f_s/(conf.ofdm.ncarrier*conf.ofdm.spacing);

conf.sc.txpulse_length = 20*conf.sc.os_factor;
conf.sc.txpulse        = rrc(conf.sc.os_factor,0.22,conf.sc.txpulse_length);

rng(0);

%% ----------------- CONSISTENCY CHECKS -----------------
bitsPerSym        = conf.modulation_order;          
N                 = conf.ofdm.ncarrier;             
bitsPerOfdmSymbol = bitsPerSym * N;                

if mod(conf.nbits, bitsPerOfdmSymbol) ~= 0
    error('conf.nbits must be a multiple of %d', bitsPerOfdmSymbol);
end

nDataOfdmSym = conf.nbits / bitsPerOfdmSymbol;
fprintf('Number of OFDM data symbols: %d\n', nDataOfdmSym);

%% ----------------- GENERATE DATA & TRANSMIT -----------------
txbits = randi([0 1], conf.nbits, 1);

[txsignal, conf] = txofdm(txbits, conf);

rawtxsignal = [ zeros(conf.f_s,1) ; txsignal ; zeros(conf.f_s,1) ];
rawtxsignal = [ rawtxsignal  zeros(size(rawtxsignal)) ];

rxsignal = channel_emulator(rawtxsignal(:,1), conf);

%% ----------------- RX BASE (task1_basic) -----------------
conf_base = conf;
conf_base.rx_mode = "task1";

[rxbits_base, conf_base] = rxofdm(rxsignal, conf_base);

if length(rxbits_base) ~= length(txbits)
    error('Length mismatch (base): rxbits=%d, txbits=%d', length(rxbits_base), length(txbits));
end

%% ----------------- RX TRACKING (task2_phase) -----------------
conf_tr = conf;
conf_tr.rx_mode = "task2";
conf_tr.ofdm.alpha = 1;

[rxbits_tr, conf_tr] = rxofdm(rxsignal, conf_tr);

if length(rxbits_tr) ~= length(txbits)
    error('Length mismatch (tracking): rxbits=%d, txbits=%d', length(rxbits_tr), length(txbits));
end

%% ----------------- BER PER OFDM SYMBOL -----------------
txbits_mat    = reshape(txbits,      bitsPerOfdmSymbol, nDataOfdmSym);
rxbits_base_m = reshape(rxbits_base, bitsPerOfdmSymbol, nDataOfdmSym);
rxbits_tr_m   = reshape(rxbits_tr,   bitsPerOfdmSymbol, nDataOfdmSym);

BER_sym_base = zeros(1, nDataOfdmSym);
BER_sym_tr   = zeros(1, nDataOfdmSym);

for k = 1:nDataOfdmSym
    BER_sym_base(k) = mean(rxbits_base_m(:,k) ~= txbits_mat(:,k));
    BER_sym_tr(k)   = mean(rxbits_tr_m(:,k)   ~= txbits_mat(:,k));
end

BER_tot_base = mean(rxbits_base ~= txbits);
BER_tot_tr   = mean(rxbits_tr   ~= txbits);

fprintf('Global BER base      : %.3e\n', BER_tot_base);
fprintf('Global BER tracking  : %.3e\n', BER_tot_tr);

%% ----------------- Cumulative BER VS N-SYMBOL -----------------
BER_cum_base = zeros(1, nDataOfdmSym);
BER_cum_tr   = zeros(1, nDataOfdmSym);

for k = 1:nDataOfdmSym
    BER_cum_base(k) = mean(rxbits_base(1:k*bitsPerOfdmSymbol) ~= txbits(1:k*bitsPerOfdmSymbol));
    BER_cum_tr(k)   = mean(rxbits_tr(1:k*bitsPerOfdmSymbol)   ~= txbits(1:k*bitsPerOfdmSymbol));
end

%% ----------------- PLOT: BER PER OFDM SYMBOL -----------------
outdir = 'plots_task2';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

fig1 = figure;
semilogy(1:nDataOfdmSym, BER_sym_base, '-o','LineWidth',1.2); hold on;
semilogy(1:nDataOfdmSym, BER_sym_tr,   '-x','LineWidth',1.2);
grid on;
xlabel('OFDM data symbol index');
ylabel('BER per symbol');
title(sprintf('Channel ID = %d, SNR = %d dB (per-symbol BER)', ...
               conf.emulator_idx, conf.emulator_snr));
legend('Task1 RX','Task2 RX','Location','best');

if save
    saveas(fig1, fullfile(outdir, sprintf('ber_ofdm_symbol_ch%d.png', conf.emulator_idx)));
end

%% ----------------- PLOT: Cumulative BER -----------------
fig2 = figure;
semilogy(1:nDataOfdmSym, BER_cum_base, '-o','LineWidth',1.2); hold on;
semilogy(1:nDataOfdmSym, BER_cum_tr,   '-x','LineWidth',1.2);
yline(0.01, '--', 1.2);
grid on;
xlabel('Number of OFDM data symbols used (from start)');
ylabel('Cumulative BER');
title(sprintf('Channel ID = %d, SNR = %d dB (cumulative BER)', ...
               conf.emulator_idx, conf.emulator_snr));
legend('Task1 RX','Task2 RX','Location','southwest');

if save
    saveas(fig2, fullfile(outdir, sprintf('cumulative_ber_ch%d.png', conf.emulator_idx)));
end

%% ----------------- Threshold -----------------
th = 0.01;

idx_good_base = find(BER_cum_base < th);
idx_good_tr   = find(BER_cum_tr   < th);

if ~isempty(idx_good_base)
    fprintf('Base RX: BER cum < %.2f up to OFDM symbol #%d\n', th, max(idx_good_base));
else
    fprintf('Base RX: BER cum < %.2f never reached\n', th);
end

if ~isempty(idx_good_tr)
    fprintf('Tracking RX: BER cum < %.2f up to OFDM symbol #%d\n', th, max(idx_good_tr));
else
    fprintf('Tracking RX: BER cum < %.2f never reached\n', th);
end

%% ----------------- RX TRACKING Alpha comparison -----------------
alpha_values = 0:0.1:1; % Define alpha values from 0 to 1 with a step of 0.1
nAlpha = length(alpha_values);
BER_cum_tr_alpha = zeros(nAlpha, nDataOfdmSym); % Initialize cumulative BER for each alpha

for a = 1:nAlpha
    conf_tr = conf;
    conf_tr.rx_mode = "task2";
    conf_tr.ofdm.alpha = alpha_values(a); % Set the current alpha value

    [rxbits_tr, conf_tr] = rxofdm(rxsignal, conf_tr);

    if length(rxbits_tr) ~= length(txbits)
        error('Length mismatch (tracking): rxbits=%d, txbits=%d', length(rxbits_tr), length(txbits));
    end

    % Calculate cumulative BER for the current alpha
    for k = 1:nDataOfdmSym
        BER_cum_tr_alpha(a, k) = mean(rxbits_tr(1:k*bitsPerOfdmSymbol) ~= txbits(1:k*bitsPerOfdmSymbol));
    end
end

%% ----------------- PLOT: Cumulative BER for different alpha values -----------------
fig3 = figure;
hold on;
for a = 1:nAlpha
    semilogy(1:nDataOfdmSym, BER_cum_tr_alpha(a, :), '-o', 'LineWidth', 1.2, 'DisplayName', sprintf('Alpha = %.1f', alpha_values(a)));
end
yline(0.01, '--', '1% BER');
set(gca, 'YScale', 'log');
grid on;
xlabel('Number of OFDM data symbols used (from start)');
ylabel('Cumulative BER');
title(sprintf('Cumulative BER vs N-Symbol for different alpha values (Channel ID = %d, SNR = %d dB)', ...
               conf.emulator_idx, conf.emulator_snr));
legend('Location', 'southeast');

if save
    saveas(fig3, fullfile(outdir, sprintf('cumulative_ber_alpha_ch%d.png', conf.emulator_idx)));
end
