function symbols = mapper_QPSK(bits)
% MAPPER_QPSK QPSK mapper (Gray) 
% Mapping:
%   [b1 b2] = [0 0] ->  1 + 1j
%   [b1 b2] = [0 1] -> -1 + 1j
%   [b1 b2] = [1 0] -> -1 - 1j
%   [b1 b2] = [1 1] ->  1 - 1j

    bits = bits(:).';

    if mod(length(bits), 2) ~= 0
        error('mapper_QPSK: length(bits) must be a multiple of 2');
    end

    % reshape in [Nsym x 2]: [b1 b2]
    bits2 = reshape(bits, 2, []).';
    b1 = bits2(:,1);
    b2 = bits2(:,2);

    symIdx = b1*2 + b2;

    % QPSK constellation
    QPSK = [1+1j, -1+1j, -1-1j, 1-1j] / sqrt(2);

    symbols = QPSK(symIdx+1).';
end
