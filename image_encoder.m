function bits = image_encoder(img)
    % image_encoder encodes a grayscale image into a binary representation.
    % 
    % Inputs:
    %   img - A grayscale image represented as a uint8 matrix.
    %          The image must be in grayscale format (single channel).
    %
    % Outputs:
    %   bits - A column vector containing the binary representation of the image pixels.
    %          Each pixel is represented as an 8-bit binary number.
    
    img = uint8(img);
    
    % Convert to column vector
    pixels = img(:);
    
    % Convert each pixel to 8-bit binary
    bits = double(int2bit(pixels.', 8));

    % Return as a single column vector
    bits = bits(:);
end
