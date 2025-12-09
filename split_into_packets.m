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
