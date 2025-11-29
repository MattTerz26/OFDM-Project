function [preamble] = preamble_generate(length)
% preamble_generate() 
% input : length: a scaler value, desired length of preamble.
% output: preamble: preamble bits
preamble = zeros(length, 1);
LFSR_state = ones(8,1);
for i = 1: length
    %get new bit 
    new = mod(sum(LFSR_state([4 5 6 8])),2);
    %output one bit
    preamble(i) = LFSR_state(end);
    %shift all bits
    LFSR_state(2:end) = LFSR_state(1:end-1);
    % input new bit
    LFSR_state(1) = new ;   
end
end