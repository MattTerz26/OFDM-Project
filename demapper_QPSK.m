function b = demapper_QPSK(symbols)
% demapper_QPSK - Demodulates QPSK symbols to binary representation.
%
% Syntax:  b = demapper_QPSK(symbols)
%
% Inputs:
%    symbols - A vector of complex QPSK symbols (1 + 1j, 1 - 1j, -1 + 1j, -1 - 1j)
%
% Outputs:
%    b - A binary vector representing the demapped bits, where:
%        [ 1 + 1j] -> [0 0]
%        [ 1 - 1j] -> [0 1]
%        [-1 + 1j] -> [1 0]
%        [-1 - 1j] -> [1 1]
%
% Example:
%    b = demapper_QPSK([1 + 1j, -1 - 1j]);

    sym = symbols(:);       % Ensure symbols are in column format
    
    % Determine the binary values based on the real and imaginary parts
    b0 = real(sym) < 0;     % First bit based on the real part
    b1 = imag(sym) < 0;     % Second bit based on the imaginary part
    
    B = [b0 b1].';          % Create a 2 x Ns matrix
    b = B(:);               % Reshape to (2*Ns) x 1 vector
end