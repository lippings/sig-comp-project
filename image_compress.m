function [] = image_compress(input_file, output_file, method, plot_images)
% Compress image and display required information as instructed by the
% project work instructions.
% method has the following options:
%   'GR' (dafualt) - Golomb-Rice (no blocking)
%   'GRBlock' - Golomb-Rice with blocking
%   'Huffman' - Huffman
if nargin < 3
    method = 'GR';
end

if nargin < 4
     plot_images = true;
end

% Read image
A = imread(input_file);
A = double(A);

if plot_images
    % Plot histogram
    h = histogram(A(:), 'Normalization', 'probability');
    ylabel('Frequency');
    xlabel('Pixel intensity');

    % Calculate entropy
    probs = h.Values(h.Values ~= 0);
    entropy = -sum(probs.*log2(probs))
end

% Crop center and save
B = A(50:249, 50:249);
imname = ['My' input_file];
imwrite(uint8(B), imname);

if plot_images
    % Plot image with cropped rectangle
    figure();
    imshow(uint8(A));
    rectangle('Position', [50, 50, 200, 200], 'EdgeColor', 'r', 'LineWidth', 2);
    
    disp([' Reference codelength: ' num2str(8*numel(B))]);
end

e = decorrelate(B, 'column', plot_images);
use_char = false; % For testing dec2bin (binary as char array) vs de2bi (binary as int array)
if strcmp(method, 'GRBlock')
    disp('Golomb-Rice with blocking');
    L = [50:50:2000];
    CL = zeros(size(L));
    for i = 1:length(L)
        [~, cl] = optimal_p(e, L(i));
        CL(i) = cl;
    end
    [~, ind] = min(CL);
    block_size = L(ind);
    bpp_CL = CL / numel(B)
    bpp = bpp_CL;
    
    disp(['Optimal block size: ' num2str(block_size)]);
    
    tic;
    bitstream = compress_blockwise(e, block_size);
    write_binary(output_file, bitstream, use_char);
    compress_time = toc;
    clear bitstream;
    
    tic;
    bitstream = read_binary(output_file, use_char);
    recon_e = decompress_blockwise(bitstream, block_size);
    decompress_time = toc;
    
    codelen = length(bitstream);
elseif strcmp(method, 'GR')
    disp('Golomb-Rice without blocking');
    
    [~, cl] = optimal_p(e);
    bpp = cl / numel(B);
    
    tic;
    bitstream = compress_blockwise(e, nan); % nan for no blocking
    write_binary(output_file, bitstream, use_char);
    compress_time = toc;
    clear bitstream;

    tic;
    bitstream = read_binary(output_file, use_char);
    recon_e = decompress_blockwise(bitstream, nan); % nan for no blocking
    decompress_time = toc;
    
    codelen = length(bitstream);
elseif strcmp(method, 'Huffman')
    disp('Huffman coding');
    [syms, probs] = get_symbol_probs(e);
    
    ent = -sum(probs.*log2(probs));
    ideal_cl = numel(e)*ent;
    disp(['Ideal codelength (with BPP=entropy=' num2str(ent) '): ' num2str(ideal_cl)]);
    
    tic;
    [dict, bpp] = huffmandict(syms, probs);
    bitstream = huffmanenco(e, dict);
    write_binary(output_file, bitstream, use_char);
    compress_time = toc;
    clear bitstream;
    
    tic;
    bitstream = read_binary(output_file, use_char);
    recon_e = huffmandeco(bitstream, dict);
    recon_e = recon_e(1:end-1); % For some reason, an extra number appears at the end
    decompress_time = toc;
    
    n = length(syms); % Number of leaves
    codelen = length(bitstream) + 2*n - 1; % 2n-1 bits for encoding a binary tree
else
    disp(['Unrecognized method name ' method]);
    return
end

comp_filesize = human_readable_file_size(output_file);

disp(['Compression time: ' num2str(compress_time) 's Decompression time: ' num2str(decompress_time) 's']);
disp(['Compressed size: ' comp_filesize ' (Codelength ' num2str(codelen) ')']);

if length(bpp) == 1
    disp(['Bits per pixel: ' num2str(bpp)]);
else
    figure();
    plot(L, bpp);
    ylabel('Bits per pixel');
    xlabel('Block size');
%     title('BBP by block size');
end

disp(['Total compression error (should be 0): ' num2str(sum(e-recon_e, 'all'))]);
end


function [symbols, probabilities] = get_symbol_probs(values)
N = length(values);

symbols = unique(values);
counts = zeros(size(symbols));

for i = 1:length(symbols)
    counts(i) = sum(values == symbols(i));
end

probabilities = counts / N;
end


function [] = write_binary(filename, bitstream, convert_to_int)
% Write binary file. filename must have the extension .bin
if nargin < 3
    convert_to_int = false;
end

if convert_to_int
    bitstream = bitstream - '0';
end

fid = fopen(filename, 'w');
fwrite(fid, bitstream, 'ubit1');
fclose(fid);
end


function [bitstream] = read_binary(filename, convert_to_char)
% Read binary file. filename must have the extension .bin
if nargin < 2
    convert_to_char = false;
end
    
fid = fopen(filename, 'r');
bitstream = fread(fid, inf, 'ubit1');
fclose(fid);

if convert_to_char
    bitstream = num2str(bitstream)';
else
    bitstream = bitstream';
end
end


function [values] = decompress_blockwise(bitstream, block_size)
bs_ind = 1;

bs_nind = bs_ind+32-1;
array_length = binary2dec(bitstream(bs_ind:bs_nind));
bs_ind = bs_nind+1;

if ~isnan(block_size)
    bs_nind = bs_ind+16-1;
    block_size = binary2dec(bitstream(bs_ind:bs_nind));
    bs_ind = bs_nind+1;
end

values = zeros(1, array_length);

if ~isnan(block_size)
    for val_ind = 1:block_size:array_length
        if val_ind+block_size > array_length
            val_nind = array_length;
        else
            val_nind = val_ind + block_size - 1;
        end

        [block_values, bs_ind] = decompress_values(bitstream, bs_ind, val_nind-val_ind+1);
        
        values(val_ind:val_nind) = block_values;
    end
else
    values = decompress_values(bitstream, bs_ind, array_length);
end
end


function [values, end_ind] = decompress_values(bitstream, start_ind, array_length)
% Decompress bitstream into an array of values
ind = start_ind;

% p value
nind = ind+3-1;
p = binary2dec(bitstream(ind:nind));
ind = nind+1;
base = 2^p;

signs = zeros(array_length, 1);
quotients = zeros(array_length, 1);
if p == 0
    remainders = zeros(array_length, 1);
    for i = 1:array_length
        nind = ind;
        sign = binary2dec(bitstream(ind:nind));
        ind = nind+1;

        [~, quotient, nind] = read_unary(bitstream, ind);
        ind = nind+1;
        
        signs(i) = sign;
        quotients(i) = quotient;
    end
else
    remainders_bin = zeros(array_length, p);
    for i = 1:array_length
        nind = ind;
        sign = bitstream(ind:nind);
        ind = nind+1;
        
        nind = ind+p-1;
        remainder_bin = bitstream(ind:nind);
        ind = nind+1;
        
        [~, quotient, nind] = read_unary(bitstream, ind);
        ind = nind+1;
        
        signs(i) = sign;
        remainders_bin(i, :) = remainder_bin;
        quotients(i) = quotient;
    end
    remainders = binary2dec(remainders_bin);
end
signs(signs == 1) = -1;
signs(signs == 0) = 1;

values = signs.*(quotients*base+remainders);
values = values';
end_ind = ind;
end


function [bitstream] = compress_blockwise(values, block_size)
[p, codelength] = optimal_p(values, block_size);

array_length = length(values);

if ~isnan(block_size)
    bitstream = num2str(zeros(codelength + 32 + 16, 1))'; % 32 bits for array_length, 16 for block_size
else
    bitstream = num2str(zeros(codelength + 32, 1))';
end

ind = 1;

nind = ind+32-1;
bitstream(ind:nind) = dec2binary(array_length, 32);
ind = nind+1;

if ~isnan(block_size)
    nind = ind+16-1;
    bitstream(ind:nind) = dec2binary(block_size, 16);
    ind = nind+1;
end

idx = 1;
if ~isnan(block_size)
    for i = 1:block_size:array_length
        if i+block_size > array_length
            block_values = values(i:end);
        else
            block_values = values(i:i+block_size-1);
        end
        
        compressed_block = compress_array(block_values, p(idx));
        idx = idx + 1;

        nind = ind+length(compressed_block)-1;
        bitstream(ind:nind) = compressed_block;
        ind = nind+1;
    end
else
    bitstream(ind:end) = compress_array(values, p);
end
end


function [bitstream] = compress_array(values, p)
% Compress array of values using GR encoding
if nargin < 2
    [p, codelength] = optimal_p(values);
else
    codelength = GR_estimation(values, p);
end

array_length = length(values);

bitstream = zeros(1, codelength);

ind = 1;

% Add p value to the beginning of bitstream
nind = ind+3-1; % Next index
bitstream(ind:nind) = dec2binary(p, 3);
ind = nind+1;

% Do as much as possible with matrix operations.
signs = zeros(size(values));
signs(values < 0) = 1;

base = 2^p;
if p == 0
    quotients = abs(values) / base;
    
    % Write rest of bitstream in loop
    for i = 1:array_length
        quotient = quotients(i);
        quotient_unary = dec2unary(quotient);
        sign_bin = signs(i);
        
        nind = ind+1+quotient;
        bitstream(ind:nind) = [sign_bin quotient_unary];
        ind = nind+1;
    end
else
    remainders = mod(abs(values), base);
    remainders_bin = dec2binary(remainders, p);
    quotients = (abs(values)-remainders) / base;

    % Write rest of bitstream in loop
    for i = 1:array_length
        quotient = quotients(i);
        quotient_unary = dec2unary(quotient);
        remainder_bin = remainders_bin(i, :);
        sign_bin = signs(i);

        nind = ind+1+p+quotient;
        bitstream(ind:nind) = [sign_bin remainder_bin quotient_unary];
        ind = nind+1;
    end
end

end


function [e] = decorrelate(img, unroll_method, plot_images)
% Decorrelate image using MAP
if nargin < 2
    unroll_method = 'column';
end

if nargin < 3
    plot_images = false;
end
E = zeros(size(img));

[h, w] = size(img);
for i = 1:h
    for j = 1:w
        z = img(i, j);
        if j == 1
            z_hat = 0;
        else
            n = nan;
            w = nan;
            nw = nan;
            
            if i > 1
                n = img(i-1, j);
            end
            
            if j > 1
                w = img(i, j-1);
            end
            
            if i > 1 && j > 1
                nw = img(i-1, j-1);
            end
            
            vals = [n w n+w-nw];
            
            z_hat = round(median(vals(~isnan(vals))));
        end
        
        E(i, j) = z - z_hat;
    end
end

if plot_images
    % Show result
    figure();
    imagesc(E);
end

if strcmp(unroll_method, 'zig-zag')
    e = zigzag(E);
elseif strcmp(unroll_method, 'row')
    E = E';
    e = E(:);
else
    e = E(:)';
end
end


function [pstar, best_cl] = optimal_p(values, block_size)
% Estimate the optimal value p for Goulomb-Rice encoding for the given
% values.
% Input values are expected to range between -2^8, 2^8
if nargin < 2
    block_size = nan;
end

array_length = length(values);

if isnan(block_size)
    code_lengths = zeros(1, 8);
    p_cand = 0:7;
    for p = p_cand
        cl = GR_estimation(values, p);
        code_lengths(p+1) = cl;
    end
    
    [best_cl, ind] = min(code_lengths);
    pstar = p_cand(ind);
else
    pstar = zeros(1, ceil(array_length/block_size));
    best_cl = 0;
    idx = 1;
    for i = 1:block_size:array_length
        if i+block_size > array_length
            block_values = values(i:end);
        else
            block_values = values(i:i+block_size-1);
        end
        
        code_lengths = zeros(1, 8);
        p_cand = 0:7;
        for p = p_cand
            cl = GR_estimation(block_values, p);
            code_lengths(p+1) = cl;
        end
        
        [best_cl_block, ind] = min(code_lengths);
        pstar_block = p_cand(ind);
        
        pstar(idx) = pstar_block;
        idx = idx + 1;
        best_cl = best_cl + best_cl_block;
    end
end
end


function [codelength] = GR_estimation(values, p)
% Evaluate code length using Goulomb Rice for the given values and the
% given p value. Also includes the 3 bits needed for the p value (but no
% other meta information).
base = 2^p;
codelength = sum(1 + p + floor(abs(values) / base) + 1) + 3;
end


function [unary_value, dec_value, end_ind] = read_unary(bitstream, start_ind)
% Read the next unary number from bitstream starting at start_ind
end_ind = find_next_1(bitstream, start_ind);
unary_value = bitstream(start_ind:end_ind);
dec_value = length(unary_value) - 1;
end


function [end_ind] = find_next_1(bitstream, start_ind)
brk = [false true];

end_ind = start_ind;
while ~brk(bitstream(end_ind)+1)
% while ~brk((bitstream(end_ind)-'0')+1)
    end_ind = end_ind + 1;
end
end


function [dec_value] = binary2dec(binary_value)
% dec_value = bin2dec(binary_value);
dec_value = bi2de(binary_value, 'left-msb');
end


function [binary_value] = dec2binary(dec_value, nb_bits)
% binary_value = dec2bin(dec_value, nb_bits);
binary_value = de2bi(dec_value, nb_bits, 'left-msb');
end


function [unary_value] = dec2unary(dec_value)
% Convert decimal value to unary value
% unary_value = [repmat('0', 1, dec_value), '1'];
unary_value = [zeros(1, dec_value) 1];
end

