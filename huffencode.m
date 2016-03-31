function [zipped, info] = huffencode(vector)
% ������������ uint8 ��ʽ
% info ���ؽ�����Ҫ�Ľṹ��Ϣ
% info.pad ����ӵı�����
% info.huffcodes �� Huffman ����
% info.rows ��ԭʼͼ������
% info.cols ��ԭʼͼ������
% info.length ��ԭʼͼ�����ݳ���
% info.maxcodelen ������볤
% info.totalCodeNo ���ܵ�����

% if ~isa(vector, 'uint8')
%     error('input argument must be a uint8 vector');
% end
[m, n] = size(vector);
vector=double(vector(:));
tic;
[totalNum,tempNum]=checkNum(vector);        %   totalNum��vector�еĲ�ͬ���ֵĸ���
eleTable=0;
for i=1:totalNum
    eleTable(i)=tempNum{i}(1);              %   ��ȡtempNum�е�ÿ��Ԫ���ĵ�һ��Ԫ�أ��ŵ�eleTable��
end


vector = vector(:)';
f = frequency(vector);      %��������ų��ֵĸ���
symbols = find(f~=0);
f = f(symbols);
[f, sortindex] = sort(f);    %�����Ű��ճ��ֵĸ��ʴ�С����
symbols = symbols(sortindex);
len = length(symbols);
symbols_index = num2cell(1:len);
codeword_tmp = cell(len, 1);


%=========================================================================
%�������� generate the codeword as the 52 bits of a double
%   simbols_index ��Ƶ��f��������29��30line�����µĽڵ㣨����Ҫע�⣺1.�µĽڵ��ж�Ӧ�Ÿ���������������ֵ�ļ���2.f��֮ǰ
%   ��10 line�����Ѿ��Ź���ģ���������simbols_index��ʼ��ʱ����1�� len��19 line����������ǰ���������޹ء���
%   codeword_temp�洢���֣��洢ʱ����12����len��˳��洢�ģ���Ϊ��ֵ���Ϊ��codeword_temp(����)
%=========================================================================
% ���� Huffman �����õ����ֱ����
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

for index = 1:length(vector)    %�õ�����ͼ�����б�����
  %  (find(vector(index)==eleTable);
    len = len + length(codeword{find(vector(index)==eleTable)});  
end
string = repmat(uint8(0), 1, len);
pointer = 1;
for index = 1:length(vector)       %������ͼ����б���
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
    
%����ѹ��������
cols = length(string)/8;
string = reshape(string, 8, cols);
weights = 2.^(0: 7);
zipped = uint8(weights * double(string));
    
%���洢��һ��ϡ�����
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

%��д����ʱ����Ľṹ��Ϣ
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
    
%����addnode��ӽڵ�
function codeword_new = addnode(codeword_old, item)
codeword_new = cell(size(codeword_old));
for index = 1:length(codeword_old)
    codeword_new{index} = [item codeword_old{index}];
end

%����frequency��������ų��ֵĸ���
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
f = f./len;   %��һ��
    


    