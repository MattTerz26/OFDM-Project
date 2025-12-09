function symbols = mapper_QPSK(bits)
% MAPPER_QPSK QPSK mapper
% Mapping:
%   [b0 b1] = [0 0] ->  1 + 1j
%   [b0 b1] = [0 1] ->  1 - 1j
%   [b0 b1] = [1 0] -> -1 + 1j
%   [b0 b1] = [1 1] -> -1 - 1j

    % bits = bits(:).';
    % 
    % if mod(length(bits), 2) ~= 0
    %     error('mapper_QPSK: length(bits) must be a multiple of 2');
    % end
    % 
    % % reshape in [Nsym x 2]: [b1 b2]
    % bits2 = reshape(bits, 2, []).';
    % b0 = bits2(:,1);
    % b1 = bits2(:,2);
    % 
    % symIdx = b0*2 + b1;
    % 
    % % QPSK constellation
    % QPSK = [1+1j, 1-1j, -1+1j, -1-1j] / sqrt(2);
    % 
    % symbols = QPSK(symIdx+1).';

% QPSK mapper for bit matrices and vectors: each column = one OFDM symbol
% Input:  bits (2N x Ns) of {0,1}
% Output: sym  (N x Ns) complex matrix

    % Size
    [rows, ~] = size(bits);
    if mod(rows, 2) ~= 0
        error("QPSK mapper: number of rows must be even (2 bits per symbol).");
    end
    
    % Split bit planes
    b0 = bits(1:2:end, :);   % MSB
    b1 = bits(2:2:end, :);   % LSB

    % Gray QPSK mapping
    % 0 -> +1, 1 -> –1

    I = 1 - 2*b0;
    Q = 1 - 2*b1;

    symbols = ( double(I) + 1j.*double(Q)) / sqrt(2);   % normalize to unit power

end
