function rxsignal = channel_known(txsignal, conf)
% CHANNEL_EMULATOR  Canale noto per debug della pipeline OFDM
%
%   rxsignal = channel_emulator(txsignal, conf)
%
%   txsignal : segnale reale passband a frequenza conf.f_s
%   conf    : struttura di configurazione
%
%   Usa conf.known_channel_case per scegliere il canale
%   Usa conf.emulator_snr (in dB) per il rumore AWGN

    txsignal = txsignal(:);             % colonna

    % Scegli il tipo di canale noto
    if ~isfield(conf, "known_channel_case")
        conf.known_channel_case = 1;
    end

    Fs = conf.f_s;                      % solo per chiarezza, qui non serve molto

    switch conf.known_channel_case
        case 1
            % Caso 1: canale piatto, solo rumore
            h = 1;

        case 2
            % Caso 2: canale a 3 cammini, tutti dentro il CP
            % ritardi in campioni audio
            d = [0 3 7];                % posizioni dei tap
            g = [1 0.5 0.3];           % ampiezze dei tap
            L = d(end) + 1;
            h = zeros(L,1);
            h(d+1) = g;

        case 3
            % Caso 3: 5 tap con decadimento esponenziale
            L = 8;
            n = (0:L-1).';
            h = exp(-n/2);              % decadimento
            h = h ./ norm(h);           % normalizza un minimo

        otherwise
            error("channel_emulator: known_channel_case sconosciuto");
    end

    % Applica il canale
    rx_clean = conv(txsignal, h);

    % Porta la lunghezza come in ingresso
    rx_clean = rx_clean(1:length(txsignal));

    % Aggiungi AWGN in funzione di conf.emulator_snr
    if isfield(conf, "emulator_snr") && isfinite(conf.emulator_snr)
        SNRdB  = conf.emulator_snr;
        SNRlin = 10^(SNRdB/10);

        sigPow   = mean(rx_clean.^2);
        noisePow = sigPow / SNRlin;

        noise = sqrt(noisePow) * randn(size(rx_clean));  % rumore reale
        rxsignal = rx_clean + noise;
    else
        rxsignal = rx_clean;
    end
end
