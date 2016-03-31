function [ sig, head ] = daeread( inputfile )
%DAEREAD Summary of this function goes here
%   Detailed explanation goes here

%% read file
% open file
fid = fopen(inputfile, 'rb');
% file tags
head.stag = fread(fid, 1, 'ubit8');% HEX:0x91
head.suffix = fread(fid, 3, 'ubit8');% DAE, HEX:0x44 0x41 0x45
head.imgDims = fread(fid, 1, 'ubit8');% image dims
head.imgSize = fread(fid, head.imgDims, 'ubit16');
head.comDims = fread(fid, 1, 'ubit8');% compressed image dims
head.comSize = fread(fid, head.comDims, 'ubit16');
head.numSymbols = fread(fid, 1, 'ubit32');
% dict
symbols = zeros(head.numSymbols,1);
for i = 1:head.numSymbols
    symbols(i) = fread(fid, 1, 'ubit8');
end
bitcount = fread(fid, 1, 'ubit8');
counts = zeros(head.numSymbols,1);
if bitcount == 16
    counts = fread(fid, head.numSymbols, 'ubit16');
elseif bitcount == 32
    counts = fread(fid, head.numSymbols, 'ubit32');
end
bslen = fread(fid, 1, 'ubit64');
sig_encoded = fread(fid, bslen, 'ubit1');

fclose(fid);

prob = counts./sum(counts);
[dict,avglen] = huffmandict(symbols,prob);
sig = huffmandeco(sig_encoded, dict);
flag = 0;