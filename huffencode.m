function [zipped, info] = huffencode(vector)
% 输入和输出都是 uint8 格式
% info 返回解码需要的结构信息
% info.pad 是添加的比特数
% info.huffcodes 是 Huffman 码字
% info.rows 是原始图像行数
% info.cols 是原始图像列数
% info.length 是原始图像数据长度
% info.maxcodelen 是最大码长
% info.totalCodeNo 是总的码数

% if ~isa(vector, 'uint8')
%     error('input argument must be a uint8 vector');
% end
[m, n] = size(vector);
vector=double(vector(:));
tic;
[totalNum,tempNum]=checkNum(vector);        %   totalNum是vector中的不同数字的个数
eleTable=0;
for i=1:totalNum
    eleTable(i)=tempNum{i}(1);              %   抽取tempNum中的每个元胞的第一个元素，放到eleTable中
end


vector = vector(:)';
f = frequency(vector);      %计算各符号出现的概率
symbols = find(f~=0);
f = f(symbols);
[f, sortindex] = sort(f);    %将符号按照出现的概率大小排列
symbols = symbols(sortindex);
len = length(symbols);
symbols_index = num2cell(1:len);
codeword_tmp = cell(len, 1);


%=========================================================================
%产生码字 generate the codeword as the 52 bits of a double
%   simbols_index 是频率f重新排序（29，30line）后新的节点（这里要注意：1.新的节点中对应着该子树下子树索引值的集合2.f在之前
%   （10 line）就已经排过序的，而且这里simbols_index初始化时就是1到 len（19 line），早与以前排序索引无关。）
%   codeword_temp存储码字，存储时按照12……len的顺序存储的，因为赋值语句为：codeword_temp(索引)
%=========================================================================
% 生成 Huffman 树，得到码字编码表
while length(f)>1
    index1 = symbols_index{1};
    index2 = symbols_index{2};
    codeword_tmp(index1) = addnode(codeword_tmp(index1), uint8(0));
    codeword_tmp(index2) = addnode(codeword_tmp(index2), uint8(1));    
    f = [sum(f(1:2)),f(3:end)];
    symbols_index = [{[index1, index2]},symbols_index(3:end)];    
    [f, sortindex] = sort(f); 
    symbols_index = symbols_index(sortindex);
end
codeword = cell(totalNum, 1);
codeword(symbols) = codeword_tmp;
len = 0;

for index = 1:length(vector)    %得到整个图像所有比特数
  %  (find(vector(index)==eleTable);
    len = len + length(codeword{find(vector(index)==eleTable)});  
end
string = repmat(uint8(0), 1, len);
pointer = 1;
for index = 1:length(vector)       %对输入图像进行编码
    code = codeword{find(vector(index)==eleTable)};           
    len = length(code);
    string(pointer + (0:len-1))=code;
    pointer = pointer + len;
end
len = length(string);
pad = 8-mod(len, 8);
if pad > 0
    string = [string uint8(zeros(1, pad))];
end

info.totalCodeNo = length(string);
info.string = string;

codeword = codeword(symbols);
symLen = length(codeword);
codelen = zeros(size(codeword));
weights = 2.^(0:23);
maxcodelen = 0;
for index = 1:length(codeword)
    len = length(codeword{index});
    if len > maxcodelen;
        maxcodelen = len;
    end
    if len > 0
        code = sum(weights(codeword{index} == 1));
        code = bitset(code, len + 1);
        codeword{index} = code;
        codelen(index) = len;
    end
end
codeword = [codeword{:}];
    
%计算压缩的向量
cols = length(string)/8;
string = reshape(string, 8, cols);
weights = 2.^(0: 7);
zipped = uint8(weights * double(string));
    
%码表存储到一个稀疏矩阵
huffcodes = sparse(1, 1);
for index = 1:nnz(codeword)   % length(codeword)  %numel(codeword)
%     huffcodes(codeword(index), 1) = symbols(index);
     huffcodes(codeword(index), 1) = eleTable(symbols(index));
end
huffTIndex = find(huffcodes);
temphuffcodes = full(huffcodes);
huffT = temphuffcodes(huffTIndex);
huffT = huffT - 1;
for i = 1:length(huffTIndex)
     temp = dec2bin(huffTIndex(i));
     temp = temp(2:end);
     huffTIndex(i) = bin2dec(temp);
end

%填写解码时所需的结构信息
info.pad = pad;
info.huffcodes = huffcodes;
info.huffTIndex = huffTIndex;
info.huffT = huffT;
info.symLen = symLen;
info.ratio = cols./length(vector);
info.length = length(vector);
info.maxcodelen = maxcodelen;
info.rows = m;
info.cols = n;
    
%函数addnode添加节点
function codeword_new = addnode(codeword_old, item)
codeword_new = cell(size(codeword_old));
for index = 1:length(codeword_old)
    codeword_new{index} = [item codeword_old{index}];
end

%函数frequency计算各符号出现的概率
function f = frequency(vector)
% if ~isa(vector, 'uint8')
%     error('input argument must be a uint8 vector');
% end
[totalNum,tempNum]=checkNum(vector);
f = repmat(0, 1, totalNum);
len = length(vector);
for i=1:totalNum
    f(i)=tempNum{i}(2);
end
% for index = 0:255
%     f(index+1) = sum(vector == uint8(index));
% end
f = f./len;   %归一化
    


    