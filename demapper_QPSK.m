function b = demapper_QPSK(symbols)
% QPSK = [1+1j, -1+1j, -1-1j, 1-1j] / sqrt(2)
% bits → symbol:
%   00 ->  1+1j
%   01 -> -1+1j
%   10 -> -1-1j
%   11 ->  1-1j

    sym = symbols(:);       % col

    sR = real(sym) < 0;
    sI = imag(sym) < 0;

    bit1 = sI;
    bit2 = xor(sR, sI);

    B = [bit1 bit2].';      % 2 x Ns
    b = B(:);               % (2*Ns) x 1
end
