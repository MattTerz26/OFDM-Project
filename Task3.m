close all; clear all; clc;

save = 1;    % =1 save pictures
rng(0);


%% BASE CONFIGURATION

conf.audiosystem    = 'emulator';
conf.emulator_idx   = 2;        
conf.emulator_snr   = 25;

conf.f_c     = 8000;
conf.Nsym    = 100;
NsymTot      = conf.Nsym + 1;   % + training symbol
conf.nbits   = 512*2*conf.Nsym;   

% SC preamble
conf.sc.f_sym = 1000;          
conf.sc.nsyms = 500;     

% OFDM
conf.ofdm.bandwidth = 2000;    
conf.ofdm.ncarrier  = 512;     
conf.ofdm.cplen     = 256;     
conf.modulation_order = 2;     % QPSK (2 bits per sym)

% Audio FS
conf.f_s    = 48000;           
conf.bitsps = 16;

% Derived
conf.ofdm.spacing   = conf.ofdm.bandwidth/conf.ofdm.ncarrier;
conf.sc.os_factor   = conf.f_s/conf.sc.f_sym;
conf.ofdm.os_factor = conf.f_s/(conf.ofdm.ncarrier*conf.ofdm.spacing);

conf.sc.txpulse_length = 20*conf.sc.os_factor;
conf.sc.txpulse        = rrc(conf.sc.os_factor, 0.22, conf.sc.txpulse_length);

%% CONSISTENCY CHECKS
bitsPerSym        = conf.modulation_order;          
N                 = conf.ofdm.ncarrier;             
bitsPerOfdmSymbol = bitsPerSym * N;

if mod(conf.nbits, bitsPerOfdmSymbol) ~= 0
    error('conf.nbits must be a multiple of %d', bitsPerOfdmSymbol);
end

nDataOfdmSym = conf.nbits / bitsPerOfdmSymbol;
nOfdmSymTot  = nDataOfdmSym + 1;  % + training symbol

%% TRANSMISSION

txbits = preamble_generate(conf.nbits);     
[txsignal, conf] = txofdm(txbits, conf);

% Padding
rawtxsignal = [ zeros(conf.f_s,1) ; txsignal ; zeros(conf.f_s,1) ];
rawtxsignal = [ rawtxsignal  zeros(size(rawtxsignal)) ];



%% LOOP OVER CHANNELS

channelIDs = 2:5;
outdir = "plots_task3";
if ~exist(outdir,"dir")
    mkdir(outdir);
end

SNRdB      = 100;   % high SNR to observe channel cleanly
Nfft_delay = N;     
Nfft_freq  = N;     


for chID = channelIDs

    fprintf("\n================ CHANNEL %d ================\n", chID);

    conf.emulator_idx = chID;
    conf.emulator_snr = SNRdB;

    rng(0);
    rxsignal = channel_emulator(rawtxsignal(:,1), conf);

    
    %% RX FRONTEND (NO FRAME SYNC, PERFECT REMOVAL)

    % Remove padding
    rx_trim = rxsignal(conf.f_s+1 : end-conf.f_s);

    % Downconvert
    n = (0:length(rx_trim)-1).';
    carrier_rx = exp(-1j*2*pi*conf.f_c*n/conf.f_s);
    bb_rx = rx_trim .* carrier_rx;

    % Low-pass
    bb_rx_filt = ofdmlowpass(bb_rx, conf, conf.ofdm.bandwidth);

    % Remove preamble exactly
    Lp = length(conf.sc.txpulse);
    preamble_len = conf.sc.nsyms*conf.sc.os_factor + Lp - 1;

    ofdm_bb_os_rx = bb_rx_filt(preamble_len+1:end);

    % Resample to OFDM Fs
    ofdm_bb_rx = ofdm_rx_resampling(ofdm_bb_os_rx, conf);

    total_ofdm_len = nOfdmSymTot * (N + conf.ofdm.cplen);
    ofdm_bb_rx = ofdm_bb_rx(1:total_ofdm_len);

    % Parallelize
    rxBlocks = reshape(ofdm_bb_rx, N + conf.ofdm.cplen, nOfdmSymTot);

    % Remove CP
    rxNoCP = rxBlocks(conf.ofdm.cplen+1:end, :);

    % FFT
    Z = fft(rxNoCP, N, 1);

    
    %% CHANNEL ESTIMATION

    H_est = Z ./ conf.ofdm.txSym;
    
    %% PDP  
    h_est = ifft(H_est, Nfft_delay, 1);
    m1 = 10;      % m1 is a symbol

    pdp       = abs(h_est(:, m1)).^2;
    pdp_mean  = mean(abs(h_est).^2, 2);
    pdp       = pdp / max(pdp);
    pdp_mean  = pdp_mean / max(pdp_mean);
    pdp_dB      = 10*log10(pdp + eps);
    pdp_mean_dB = 10*log10(pdp_mean + eps);

    delayAxis = 0:Nfft_delay-1;
    thr_dB = -20;
    idx_sig = find(pdp_mean_dB > thr_dB);
    L_eff  = idx_sig(end);

    fprintf("Effective Length L_eff = %d taps\n", L_eff);
    
    fig1 = figure; hold on;

    plot(pdp_dB); 
    plot(pdp_mean_dB,'LineWidth', 2);
    yline(thr_dB,'--');
    title(sprintf("PDP – Channel %d", chID));
    xlabel("Tap index"); ylabel("Power [dB]");
    legend("Symbol 10","Mean PDP");
    grid on;

    %% FREQUENCY RESPONSE
    H_plot = fftshift(H_est,1);
    mag_sym  = abs(H_plot(:,m1));
    mag_mean = mean(abs(H_plot),2);
    ref = max([mag_sym; mag_mean]);
    H_sym_dB  = 20*log10(mag_sym/ref + eps);
    H_mean_dB = 20*log10(mag_mean/ref + eps);

    fAxis = linspace(-0.5,0.5,N);

    fig2 = figure; hold on;
    plot(fAxis,H_sym_dB);
    plot(fAxis,H_mean_dB,'LineWidth',1.5);
    grid on;
    xlabel("Normalized freq"); ylabel("|H| [dB]");
    title(sprintf("Frequency Response – Channel %d",chID));
    legend("Symbol 10","Mean");
    
    %% CHANNEL EVOLUTION
    fig3 = figure;
    imagesc(abs(H_est));
    axis xy;
    xlabel("OFDM symbol"); ylabel("Subcarrier");
    title(sprintf("Channel Evolution |H(m,n)| – Ch %d", chID));
    colorbar;

    
    %% PHASE EVOLUTION
    m0 = 50;
    fig5 = figure;
    plot(1:NsymTot, unwrap(angle(H_est(m0,:))), '-o');
    grid on;
    xlabel('OFDM symbol index n');
    ylabel('angle(H(m_0,n)) [rad]');
    title(sprintf('Phase evolution, m_0 = %d - Channel ID %d', m0, chID));

    %% SUB CARRIER CHANNEL EVOLUTION
    fig4 = figure;
    plot(1:NsymTot, abs(H_est(m0,:)), '-o');
    grid on;
    xlabel('OFDM symbol index n');
    ylabel('|H(m_0,n)|');
    title(sprintf('Time evolution in amplitude, m_0 = %d - Channel ID %d', m0, chID));
  
    %% SAVE PLOTS
    if save
        saveas(fig1, fullfile(outdir, sprintf('pdp_chID_%d.png', chID)));
        saveas(fig2, fullfile(outdir, sprintf('avg_freqresp_chID_%d.png', chID)));
        saveas(fig3, fullfile(outdir, sprintf('channel_time_ev_%d.png', chID)));
        saveas(fig4, fullfile(outdir, sprintf('channel_time_ev_detail_%d.png', chID)));
        saveas(fig5, fullfile(outdir, sprintf('phase_ev_chID_%d.png', chID)));
    end
end
