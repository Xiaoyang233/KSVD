function [ flag ] = daewrite( outfile, codes )
%DAEWRITE Writes sig to a DAE file
%   
flag = 1;

imgSize = codes.imgSize;
coefSize = codes.coefSize; 
quanMin = codes.Min;
quanMax = codes.Max;
sparsity = codes.L;
coepad = codes.coepad;
rowpad = codes.rowpad;

% ==========================================
% Huffman Table
symcoeLen = codes.symcoeLen;
symrowLen = codes.symrowLen;
maxcoeLen = codes.maxcoeLen;
maxrowLen = codes.maxrowLen;


if(maxcoeLen <= 8)
    coePercision = 'ubit8';
elseif(maxcoeLen <= 16 &&maxcoeLen > 8)
    coePercision = 'ubit16';
elseif(maxcoeLen <= 32 &&maxcoeLen > 16)
    coePercision = 'ubit32';
elseif(maxcoeLen <= 64 &&maxcoeLen > 32)
    coePercision = 'ubit64';
end

if(maxrowLen <= 8)
    rowPercision = 'ubit8';
elseif(maxrowLen <= 16 &&maxrowLen > 8)
    rowPercision = 'ubit16';
elseif(maxrowLen <= 32 &&maxrowLen > 16)
    rowPercision = 'ubit32';
elseif(maxrowLen <= 64 &&maxrowLen > 32)
    rowPercision = 'ubit64';
end

huffTcoeIndex = codes.huffTcoeIndex;
huffTcoe = codes.huffTcoe;

huffTrowIndex = codes.huffTrowIndex;
huffTrow = codes.huffTrow;
% ==========================================
% code stream
bincoeLen = codes.bincoeLen;
binrowLen = codes.binrowLen;

coebin = codes.coebin;
rowbin = codes.rowbin;
% ===========================================

%% write file
% open file
fid = fopen(outfile, 'wb');
% file tags
fwrite(fid, 145, 'ubit8');% HEX:0x91
fwrite(fid, [68 65 69], 'ubit8');% DAE, HEX:0x44 0x41 0x45
fwrite(fid, numel(imgSize), 'ubit8');% image dims
fwrite(fid, imgSize, 'ubit16');

fwrite(fid, coefSize,'ubit16');   %coeffice matrix size
% fwrite(fid, numel(comSize), 'ubit8');% compressed image dims

fwrite(fid, symcoeLen, 'ubit32');
fwrite(fid, symrowLen, 'ubit32');

fwrite(fid, bincoeLen, 'ubit32');   % coeffice bytestream byte
fwrite(fid, binrowLen, 'ubit32');   % row index bytestream byte
fwrite(fid, maxcoeLen, 'ubit8');
fwrite(fid, maxrowLen, 'ubit8');

fwrite(fid, quanMin, 'ubit8');
fwrite(fid, quanMax, 'uint16');

fwrite(fid, sparsity, 'ubit8');
fwrite(fid, coepad, 'ubit8');
fwrite(fid, rowpad, 'ubit8');
% ========================================
% Huffman Table
fwrite(fid, huffTcoeIndex, coePercision); %
fwrite(fid, huffTcoe, 'ubit8');

fwrite(fid, huffTrowIndex, rowPercision); %
fwrite(fid, huffTrow, 'ubit8');
% ========================================
% code Stream
fwrite(fid, coebin, 'ubit1');
% fwrite(fid, rowbin, 'ubit1');
% ========================================

fclose(fid);

flag = 0;

end

