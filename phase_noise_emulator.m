function [rxsignal, conf] = phase_noise_emulator(txsignal, conf)
% [PHASE_NOISE_EMULATOR] - Simula un canale con sola deriva di fase (Random Walk)

    fs = conf.f_s;
    fc = conf.f_c;
    
    % *** PARAMETRI DEL RANDOM WALK (Regolare per difficoltà) ***
    % Uso 5e-5 come varianza: abbastanza veloce da richiedere il tracking, 
    % ma non troppo veloce per un sistema funzionante.
    step_variance = 5e-3; 
    
    % Downconversion a banda base
    n = (0:length(txsignal)-1).';
    carrier_tx = exp(-1j * 2 * pi * fc * n / fs); 
    tx_bb = txsignal(:) .* carrier_tx;
    
    % --- GENERAZIONE DEL RANDOM WALK DI FASE ---
    N_samples = length(tx_bb);
    
    % 1. Genera i 'passi' del rumore di fase
    phase_steps = sqrt(step_variance) * randn(N_samples, 1); 
    
    % 2. Integra i passi per ottenere la deriva totale (Random Walk)
    total_phase_drift = cumsum(phase_steps); % Fase in radianti
    
    % 3. Applicazione della rotazione di fase
    phase_rotation = exp(1j * total_phase_drift);
    rx_bb = tx_bb .* phase_rotation;
    
    % --- Upconversion e Salvataggio per Debug ---
    carrier_rx = exp(1j * 2 * pi * fc * n / fs); 
    rxsignal = real(rx_bb .* carrier_rx);
    
    % SALVATAGGIO della VERA deriva di fase (per confronto nel plot)
    % Il segnale OFDM inizia dopo il preambolo
    len_preamble = length(conf.sc.preamble_bb);
    conf.debug.true_phase_drift = total_phase_drift(len_preamble + 1 : end); 

    % (Opzionale: aggiungere rumore bianco Gaussiano se necessario)
    % SNR_linear = 10^(conf.emulator_snr / 10);
    % ...
end