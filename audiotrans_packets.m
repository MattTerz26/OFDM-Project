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
        rawtx = [zeros(conf.f_s,1); txsignal; zeros(conf.f_s,1)];
        rawtx = [rawtx zeros(size(rawtx))];

        % --- CHANNEL ---
        rx = channel_emulator(rawtx(:,1), conf);

        % --- RX ---
        [rx_bits, conf] = rxofdm(rx, conf);

        % Append bits
        rx_bits_total = [rx_bits_total ; rx_bits(:)];
    end
end


function packets = split_into_packets(bits, packet_bits)
    % split_into_packets Splits the input bits into packets of specified size.
    %   Inputs:
    %       bits - The input bits to be split into packets.
    %       packet_bits - The number of bits per packet.
    %   Outputs:
    %       packets - A cell array containing the split packets.

    N = length(bits);
    numPackets = ceil(N / packet_bits);
    packets = cell(numPackets, 1);

    for k = 1:numPackets
        i1 = (k-1)*packet_bits + 1;
        i2 = min(k*packet_bits, N);
        pkt = bits(i1:i2);

        % Padding per ultimo pacchetto
        if length(pkt) < packet_bits
            pkt = [pkt ; zeros(packet_bits - length(pkt),1)];
        end
        
        packets{k} = pkt;
    end
end
