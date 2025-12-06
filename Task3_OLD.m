close all; 
clear all; 
clc;

save = 0; % Set this variable to control saving of plots

%% ----------------- BASE CONFIGURATION -----------------
conf.audiosystem    = 'emulator';
conf.emulator_idx   = 2;        
conf.emulator_snr   = 25;

conf.nbits   = 512*2*50;       
conf.f_c     = 8000;
conf.Nsym    = 50;
NsymTot = conf.Nsym + 1;

% SC preamble
conf.sc.f_sym = 1000;          
conf.sc.nsyms = 500;     

% OFDM
conf.ofdm.bandwidth = 2000;    
conf.ofdm.ncarrier  = 512;     
conf.ofdm.cplen     = 256;     
conf.modulation_order = 2; 

% Audio
conf.f_s    = 48000;           
conf.bitsps = 16;

% Derived
conf.ofdm.spacing   = conf.ofdm.bandwidth/conf.ofdm.ncarrier;
conf.sc.os_factor   = conf.f_s/conf.sc.f_sym;
conf.ofdm.os_factor = conf.f_s/(conf.ofdm.ncarrier*conf.ofdm.spacing);

conf.sc.txpulse_length = 20*conf.sc.os_factor;
conf.sc.txpulse        = rrc(conf.sc.os_factor, 0.22, conf.sc.txpulse_length);

rng(0);

%% ----------------- CONSISTENCY CHECKS -----------------
bitsPerSym        = conf.modulation_order;          
N                 = conf.ofdm.ncarrier;             
bitsPerOfdmSymbol = bitsPerSym * N;                

if mod(conf.nbits, bitsPerOfdmSymbol) ~= 0
    error('conf.nbits must be a multiple of %d', bitsPerOfdmSymbol);
end

nDataOfdmSym   = conf.nbits / bitsPerOfdmSymbol;
nOfdmSymTot    = nDataOfdmSym + 1;  % +1 Trainingssymbol
fprintf('Number of OFDM data symbols: %d\n', nDataOfdmSym);

%% ----------------- GENERATE DATA & TRANSMIT ONCE -----------------
txbits = preamble_generate(conf.nbits);     % not so random bits
[txsignal, conf] = txofdm(txbits, conf);

% Padding as in Task 1 / 2
rawtxsignal = [ zeros(conf.f_s,1) ; txsignal ; zeros(conf.f_s,1) ];
rawtxsignal = [ rawtxsignal  zeros(size(rawtxsignal)) ];

%% ----------------- Looping through the channels -----------------
channelIDs = 2:5;         

% output folder for the plots
outdir = 'plots_task3';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

SNRdB      = 100;          % high SNR for characterization
Nfft_delay = N;            % FFT length for CIR
Nfft_freq  = N;            % for frequency response

for chID = channelIDs

    fprintf('\n================ CHANNEL ID %d ================\n', chID);

    % channel parameterss
    conf.emulator_idx = chID;
    conf.emulator_snr = SNRdB;
    
    % send signal through channel
    rng(0);
    rxsignal = channel_emulator(rawtxsignal(:,1), conf); 
    % % Test known channel
    % conf.known_channel_case = 1;
    % rxsignal = channel_known(rawtxsignal(:,1), conf);


    % rx frontend
    bitsPerSym   = conf.modulation_order;       % 2 for QPSK
    N            = conf.ofdm.ncarrier;          % 512
    Ncp          = conf.ofdm.cplen;             % 256
    fs_audio     = conf.f_s;                    % 48 kHz
    fc           = conf.f_c;                    % 8 kHz
    os_sc        = conf.sc.os_factor;

    % Downconversion
    n = (0:length(rxsignal)-1).';
    carrier_rx = exp(-1j*2*pi*fc*n/fs_audio);
    bb_rx = rxsignal(:) .* carrier_rx;

    % Low-pass filtering
    bb_rx_filt = ofdmlowpass(bb_rx, conf, conf.ofdm.bandwidth);

    % Frame Sync 
    beginning_of_data = frame_sync_ofdm(bb_rx_filt, conf) - 2000;
    start_ofdm_idx = beginning_of_data + 1;

    if start_ofdm_idx <= 0 || start_ofdm_idx > length(bb_rx_filt)
        error('Task3: start_ofdm_idx out of range (%d).', start_ofdm_idx);
    end

    ofdm_bb_os_rx = bb_rx_filt(start_ofdm_idx:end);
    ofdm_bb_rx = ofdm_rx_resampling(ofdm_bb_os_rx, conf);
    total_ofdm_len = nOfdmSymTot * (N + Ncp);

    if length(ofdm_bb_rx) < total_ofdm_len
        error('Task3: ofdm_bb_rx too short (got %d, need %d).', ...
               length(ofdm_bb_rx), total_ofdm_len);
    end

    ofdm_bb_rx = ofdm_bb_rx(1:total_ofdm_len);
    
    % % Remove leading/trailing zeros
    % rx_trim = rxsignal(fs_audio+1 : end-fs_audio);   % keep only txsignal part
    % 
    % % Downconversion to complex baseband
    % n = (0:length(rx_trim)-1).';
    % carrier_rx = exp(-1j*2*pi*fc*n/fs_audio);       % minus sign for RX
    % bb_rx = rx_trim(:) .* carrier_rx;
    % 
    % % Low-pass filtering (extract baseband)
    % bb_rx_filt = ofdmlowpass(bb_rx, conf, conf.ofdm.bandwidth);
    % 
    % % Remove preamble in time-domain by knowing its length
    % %    Preamble length:
    % %       - conf.sc.nsyms symbols
    % %       - upsampled by os_sc
    % %       - convolved with RRC of length Lp = 2*txpulse_length + 1
    % %       => Np = nsyms*os_sc + Lp - 1
    % 
    % Lp = length(conf.sc.txpulse);
    % preamble_len = conf.sc.nsyms*os_sc + Lp - 1;
    % 
    % if length(bb_rx_filt) <= preamble_len
    %     error('rx_task1_basic: received signal too short compared to preamble length.');
    % end
    % 
    % ofdm_bb_os_rx = bb_rx_filt(preamble_len+1:end);
    % 
    % % Resample from audio Fs to OFDM Fs
    % ofdm_bb_rx = ofdm_rx_resampling(ofdm_bb_os_rx, conf);
    % 
    % % Shorten to expected length (multiple of (N+Ncp)*nOfdmSymTot )
    % total_ofdm_len = nOfdmSymTot * (N + Ncp);
    % 
    % if length(ofdm_bb_rx) < total_ofdm_len
    %     error('rx_task1_basic: ofdm_bb_rx too short (got %d, need %d).', ...
    %            length(ofdm_bb_rx), total_ofdm_len);
    % end
    % ofdm_bb_rx = ofdm_bb_rx(1:total_ofdm_len);
    

    % parallel blocks
    rxBlocks = reshape(ofdm_bb_rx, N+Ncp, nOfdmSymTot);

    % remove cp
    rxNoCP = rxBlocks(Ncp+1:end, :);       % N x nOfdmSymTot

    % back to freq domain
    Z = fft(rxNoCP, N, 1);              

    %% ----- channel estimation whole transmitted symbols                     
    txSym = conf.ofdm.txSym;
    H_est = Z ./ txSym;           % N x NsymTot

    %% ----- 1) power delay profile from H_est -----
    h_est = ifft(H_est, Nfft_delay, 1);

    h_mean        = mean(h_est, 2);               % media sui simboli OFDM
    [~, mainTap]  = max(abs(h_mean).^2);          % indice del tap principale

    % % receneter the channel
    % h_est_centered = circshift(h_est, 1-mainTap, 1);
    % 
    % pdp_lin     = mean(abs(h_est_centered).^2,2);
    % pdp_lin     = pdp_lin / max(pdp_lin);
    % thr = 10^(-25/20);
    % idx_sig = find(pdp_lin > thr);
    % true_span = idx_sig(end) - idx_sig(1);
    % fprintf("Effective time span: %d taps\n", true_span);



    m1 = 10;
    pdp         = abs(h_est(:, m1)).^2;
    pdp_mean    = mean(abs(h_est).^2, 2);
    pdp         = pdp / max(pdp + eps);                   % normalization -> pdp is between 0..1; eps prevents division by 0
    pdp_mean    = pdp_mean / max(pdp_mean + eps);
    pdp_dB      = 10*log10(pdp + eps);
    pdp_mean_dB = 10*log10(pdp_mean + eps);

    delayAxis = 0:Nfft_delay-1;
    thr_dB  = -10;
    idx_sig = find(pdp_mean_dB > thr_dB);
    L_eff   = max(idx_sig);
    
    fprintf('Channel %d: Effective Length (L_eff): %d\n', chID, L_eff);
    
    % NEW PLOT: more understandable
    fig1 = figure;
    hold on;
    plot(delayAxis, pdp_dB, 'LineWidth', 1.2);
    plot(delayAxis, pdp_mean_dB, 'LineWidth', 1.2);
    yline(thr_dB, '--');
    grid on;
    xlabel('Tap index (sample)');
    ylabel('Normalized power [dB]');
    title(sprintf('Power Delay Profile - Channel %d', chID));
    legend(sprintf('PDP Symbol %d', m1), 'PDP Mean');
    
    if save
        saveas(fig1, fullfile(outdir, sprintf('pdp_chID_%d.png', chID)));
    end

    % OLD PLOT
    % fig1 = figure;
    % hold on;
    % stem(delayAxis, pdp_dB, 'filled');
    % stem(delayAxis, pdp_mean_dB, 'filled');
    % xlabel('Tap index (sample)');
    % ylabel('Normalized power [dB]');
    % title(sprintf('Power Delay Profile - Channel ID %d', chID));
    % grid on;
    % legend(sprintf('Sample %d', m1),'Average over all samples');
    % saveas(fig1, fullfile(outdir, sprintf('pdp_chID_%d.png', chID)));

    %% ----- 2) frequency response -----
    H_plot = fftshift(H_est, 1);
    
    mag_sym  = abs(H_plot(:, m1));
    mag_mean = mean(abs(H_plot), 2);

    ref = max([mag_sym; mag_mean]);
    
    H_sym_dB  = 20*log10(mag_sym  / ref + eps);
    H_mean_dB = 20*log10(mag_mean / ref + eps);
    
    fAxis = linspace(-0.5, 0.5, N);
    
    fig2 = figure; hold on;
    plot(fAxis, H_sym_dB,  'LineWidth', 1.0);
    plot(fAxis, H_mean_dB, 'LineWidth', 1.5);
    hold off;
    xlabel('Normalized frequency');
    ylabel('Normalized |H(f)| [dB]');
    title(sprintf('Normalized frequency response - Channel ID %d', chID));
    legend(sprintf('Symbol %d', m1), 'Mean over symbols');
    grid on;
    
    if save
        saveas(fig2, fullfile(outdir, sprintf('avg_freqresp_chID_%d.png', chID)));
    end

    % OLD PLOT
    % H_plot = fftshift(H_est); % center zero frequency 
    % fAxis = linspace(-0.5, 0.5, Nfft_freq);
    % 
    % H_dB = 20*log10(abs(H_plot) + eps); 
    % H_mean_dB = mean(H_dB, 2); % average over all symbols
    % 
    % fig2 = figure; hold on; 
    % plot(fAxis, H_dB(:, m1), 'LineWidth', 1.2); 
    % plot(fAxis, H_mean_dB, 'LineWidth', 1.2); 
    % hold off; 
    % xlabel('Normalized frequency'); 
    % ylabel('|H(f)| [dB]'); 
    % title(sprintf('Average Frequency response - Channel ID %d', chID)); 
    % legend(sprintf('Response at Symbol %d', m1), 'Average Response'); 
    % grid on; 
    % saveas(fig2, fullfile(outdir, sprintf('avg_freqresp_chID_%d.png', chID)));

    %% ----- 3) Channel evolution over time
    fig3 = figure;
    imagesc(1:NsymTot, 1:N, abs(H_est));
    axis xy;
    xlabel('OFDM symbol index (n)');
    ylabel('Subcarrier index (m)');
    title(sprintf('|H(m,n)| - Channel ID %d', chID));
    colorbar;
    if save
        saveas(fig3, fullfile(outdir, sprintf('channel_time_ev_%d.png', chID)));
    end

    % sub carrier channel evolution
    m0 = floor(N/2) + 1;   % haf band subcarrer
    
    fig4 = figure;
    subplot(2,1,1);
    plot(1:NsymTot, abs(H_est(m0,:)), '-o');
    grid on;
    xlabel('OFDM symbol index n');
    ylabel('|H(m_0,n)|');
    title(sprintf('Time evolution in amplitude, m_0 = %d - Channel ID %d', m0, chID));
    
    subplot(2,1,2);
    plot(1:NsymTot, unwrap(angle(H_est(m0,:))), '-o');
    grid on;
    xlabel('OFDM symbol index n');
    ylabel('angle(H(m_0,n)) [rad]');
    title(sprintf('Time evolution in phase, m_0 = %d - Channel ID %d', m0, chID));
    if save
        saveas(fig4, fullfile(outdir, sprintf('channel_time_ev_detail_%d.png', chID)));
    end
    
    %% ----- 4) RMS delay spread value -----
    pdp_lin = pdp / sum(pdp);
    mu_tau  = sum(delayAxis(:) .* pdp_lin);
    tau2    = sum((delayAxis(:).^2) .* pdp_lin);
    rmsDelaySpread = sqrt(tau2 - mu_tau^2);

    fprintf('Channel %d: RMS delay spread (in taps): %.2f\n', chID, rmsDelaySpread);

end
