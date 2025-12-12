function [rx_bits_total, conf] = audiotrans_packets(bits_tx, conf)
    % audiotrans_packets Transmits and receives audio packets using OFDM modulation.
    %   Inputs:
    %       bits_tx - The transmitted bits to be sent.
    %       conf - Configuration structure containing parameters for transmission.
    %   Outputs:
    %       rx_bits_total - The total received bits after processing all packets.
    %       conf - Updated configuration structure after transmission and reception.

    bitsPerOFDM = conf.ofdm.ncarrier * conf.modulation_order;  % 1024
    maxOFDMsymbols = conf.ofdm.OFDMnSyms;
    bitsPerPacket = bitsPerOFDM * maxOFDMsymbols;

    packets = split_into_packets(bits_tx, bitsPerPacket);

    rx_bits_total = [];

    for k = 1:length(packets)
        fprintf("=== Packet %d / %d ===\n", k, length(packets));

        conf.nbits = length(packets{k});

        % --- TX ---
        [txsignal, conf] = txofdm(packets{k}, conf);
        rawtxsignal = [zeros(conf.f_s,1); txsignal; zeros(conf.f_s,1)];
        rawtxsignal = [rawtxsignal zeros(size(rawtxsignal))];

        % --- CHANNEL ---
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

        % --- RX ---
        [rx_bits, conf] = rxofdm(rxsignal, conf);

        % Append bits
        rx_bits_total = [rx_bits_total ; rx_bits(:)];
    end
end