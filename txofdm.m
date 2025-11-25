function [txsignal, conf] = txofdm(txbits,conf)
    % Digital Transmitter
    %
    %   [txsignal conf] = tx(txbits,conf) implements a complete transmitter
    %   using a single carrier preamble followed by OFDM data in digital domain.
    %
    %   txbits  : Information bits
    %   conf    : Universal configuration structure
    
    M = 2^conf.modulation_order; %4 in case of QPSK...
    bitsPerSym = conf.modulation_order; %2
    N = conf.ofdm.ncarrier; %512
    Ncp = conf.ofdm.cplen; %256
    
    %in 2 bit groups
    txbits_reshaped = reshape(txbits, bitsPerSym, []).';
    symIdx = txbits_reshaped(:,1)*2+txbits_reshaped(:,2); %symbol 0..3
    
    %QPSK mapping (gray)
    QPSK = [1+1j, -1+1j, -1-1j, 1-1j] / sqrt(2); %normalized
    qpskSyms = QPSK(symIdx+1).';
    
    %Symbols in OFDM matrix -> N x Nsym
    nOfdmSym = length(qpskSyms)/N;
    ofdmSymbols = reshape(qpskSyms,N,nOfdmSym); %512 x 50
    
    txBlocks = zeros(N+Ncp,nOfdmSym);
    for i = 1:nOfdmSym
        %IFFT
        s = ifft(ofdmSymbols(:,i)); %now in time domain...
        cp = s(end-Ncp+1:end);
        %combine payload and cp
        txBlocks(:,i) = [cp; s];
    end
    
    %serialize
    txsignal = txBlocks(:);



end