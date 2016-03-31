function testIm2jpg

% quality = 1;
% a = imread('../image/sample000.tif');
% a = rgb2gray(a);
% b = im2jpeg(a, quality);

% vector=[0,6,7,9,5,0,0,7,24,9,9,9,23];
% 
% [totalNum,temp]=checkNum(vector);

img = imread('../image/sample000.tif');
img = rgb2gray(img);  
quality = 1;   % Default value for quality.

m = [16 11  10  16  24  40  51  61        % JPEG normalizing array
     12  12  14  19  26  58  60  55       % and zig-zag redordering
     14  13  16  24  40  57  69  56       % pattern.
     14  17  22  29  51  87  80  62
     18  22  37  56  68  109 103 77
     24  35  55  64  81  104 113 92
     49  64  78  87  103 121 120 101
     72  92  95  98  112 100 103 99] * quality;

order = [1 9  2  3  10 17 25 18 11 4  5  12 19 26 33  ...
        41 34 27 20 13 6  7  14 21 28 35 42 49 57 50  ...
        43 36 29 22 15 8  16 23 30 37 44 51 58 59 52  ...
        45 38 31 24 32 39 46 53 60 61 54 47 40 48 55  ...
        62 63 56 64];
x = img;

figure(1);
imshow(x,[]);

[xm, xn] = size(x);                % Get input size.
x = double(x) - 128;               % Level shift input
t = dctmtx(8);                     % Compute 8 x 8 DCT matrix

% Compute DCTs of 8x8 blocks and quantize the coefficients.
y1 = blkproc(x, [8 8], 'P1 * x * P2', t, t');
y1 = blkproc(y1, [8 8], 'round(x ./ P1)', m);

y1 = im2col(y1, [8 8], 'distinct');  % Break 8x8 blocks into columns
xb = size(y1, 2);                   % Get number of blocks
y1 = y1(order, :);                   % Reorder column elements

eob = max(x(:)) + 1;               % Create end-of-block symbol
r = zeros(numel(y1) + size(y1, 2), 1);
count = 0;
for j = 1:xb                       % Process 1 block (col) at a time
   i = max(find(y1(:, j)));         % Find last non-zero element
   if isempty(i)                   % No nonzero block values
      i = 0;
   end
   p = count + 1;
   q = p + i;
   r(p:q) = [y1(1:i, j); eob];      % Truncate trailing 0's, add EOB,
   count = count + i + 1;          % and add to output vector
end

r((count + 1):end) = [];           % Delete unusued portion of r

[frequence,tempr] = frequency(r);

dict = huffmandict(tempr, frequence);

r_encoded = huffmanenco(r,dict);

y.size      = uint16([xm xn]);
y.numblocks = uint16(xb);
y.quality   = uint16(quality * 100);
y.huffman   = r_encoded;

save encodeHuff y dict;

rev = order;                          % Compute inverse ordering
for k = 1:length(order)
   rev(k) = find(order == k);
end

m = double(y.quality) / 100 * m;      % Get encoding quality.
xb = double(y.numblocks);             % Get x blocks.
sz = double(y.size);
xn = sz(2);                           % Get x columns.
xm = sz(1);                           % Get x rows.
x = huffmandeco(r_encoded, dict);              % Huffman decode.
eob = max(x(:));                      % Get end-of-block symbol

z = zeros(64, xb);   k = 1;           % Form block columns by copying
for j = 1:xb                          % successive values from x into
   for i = 1:64                       % columns of z, while changing
      if x(k) == eob                  % to the next column whenever
         k = k + 1;   break;          % an EOB symbol is found.
      else
         z(i, j) = x(k);
         k = k + 1;
      end
   end
end

z = z(rev, :);                                 % Restore order
x = col2im(z, [8 8], [xm xn], 'distinct');     % Form matrix blocks
x = blkproc(x, [8 8], 'x .* P1', m);           % Denormalize DCT
t = dctmtx(8);                                 % Get 8 x 8 DCT matrix
x = blkproc(x, [8 8], 'P1 * x * P2', t', t);   % Compute block DCT-1
x = uint8(x + 128);                            % Level shift

figure(2);
imshow(x,[]);

function[frequence,tempr] = frequency(r)
    totalNum = length(r);
    tempr = unique(r);
    frequence =zeros(length(tempr),1);
    for n = 1:length(tempr)
        Num = length( find(r == tempr(n)) );
        frequence(n) = Num/totalNum;
    end
end

end