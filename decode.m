function decodeVlaue = decode(code,max)
% decode for uniform quantization
% code :the uniform quantizated code ,each row store the each part code
% max :the max value of data
% clc;
% clear all;
% code=[1,0,1;1,1,0;0,0,1;0,1,0;0,0,0];
max=5;
Max = 0;
[n,m]=size(code);
a=zeros(n,1);
for i=1:n
    for j=1:m
        a(i)=a(i)+code(i,j)*2^(m-j);
    end
end
for j=1:m
    Max = Max+2^(m-j);
end
decodeVlaue = (a./Max).*max