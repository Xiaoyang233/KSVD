function vector = huffdecode(zipped, info)
% 函数对输入矩阵vector进行Huffman解码，返回解压后的图像数据

if ~isa(zipped, 'uint8')
    error('input argument must be be a uint8 vector');
end

%产生0，1序列，每位占一个字节
len = length(zipped);
string = repmat(uint8(0), 1, len.*8);
bitindex = 1:8;
for index = 1:len
    string(bitindex + 8.*(index-1)) = uint8(bitget(zipped(index), bitindex));
end
string = logical(string(:)');
len = length(string);
string ((len-info.pad+1):end)=[];
len = length(string);

%开始解码
weights = 2.^(0:51);
vector = repmat(double(0), 1, info.length);
vectorindex = 1;
codeindex = 1;
code = 0;
for index = 1:len
    code = bitset(code, codeindex, string(index));
    codeindex = codeindex+1;
    byte = decode(bitset(code, codeindex), info);
    if byte ~= 0
%         vector(vectorindex) = byte-1;
%         byte计算出的是行或者列的索引，而像素值是0-255，所以，要减1
        vector(vectorindex) = byte;
        codeindex = 1;
        code = 0;
        vectorindex = vectorindex + 1;
    end
end
vector = reshape(vector, info.rows, info.cols);

%函数decode返回码字对应的符号
function byte = decode(code, info)
byte = info.huffcodes(code);
