function img = image_decoder(bits, image_size)
    % image_decoder Decodes a binary representation of an image into a matrix.
    %   img = image_decoder(bits, image_size) takes a vector of bits and the
    %   desired image size as input, and returns the decoded image as a matrix.
    %   The image_size should be a two-element vector specifying the height and
    %   width of the image.
    %
    %   Inputs:
    %       bits - A vector containing the binary representation of the image.
    %       image_size - A two-element vector [height width] specifying the
    %                    dimensions of the output image.
    %
    %   Outputs:
    %       img - A matrix representing the decoded image.
    
    H = image_size(1);
    W = image_size(2);
    
    % Number of pixels expected
    Npix = H*W;
    
    if length(bits) ~= Npix*8
        error('Wrong number of bits for image size');
    end
    
    % Reshape bits into [Npix x 8]
    bits = reshape(bits, Npix, 8);
    
    % Convert to decimal uint8
    pixels = bit2int(bits, 8);
    
    % Reshape into image
    img = uint8(reshape(pixels, H, W));
end
