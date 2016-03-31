function[sqnr,code,a_quan]=Acode(a,n)
%a=[20*randn(1,20),randn(1,480)];
% n=16;

%A率非线性
A = 87.56;

amax=max(abs(a));
c = zeros(size(a));
for i = 1:length(a);
    if(abs((a(i)/amax))<=1/A)
        c(i) = A*(a(i)/amax)/(1+log(A));
    end
    if(abs(a(i)/amax)>1/A)
        c(i) = sign(a(i))*(1+log(A*abs(a(i)/amax)))/(1+log(A));
    end
end
        
%均匀量化
c_quan=c;
b_quan=c_quan;
d=2/n;%量化间隔
q=d.*[0:n-1];
q=q-((n-1)/2)*d;%量化电平
for i=1:n
%定位第i个量化间隔码子
    c_quan(find((q(i)-d/2<=c_quan) & (c_quan<=q(i)+d/2)))=...
    q(i).*ones(1,length(find((q(i)-d/2<=c_quan)&(c_quan<=q(i)+d/2))));
%赋值为相应的量化电平
    b_quan(find(c_quan==q(i))) =(i-1) .* ones(1,length(find(c_quan==q(i))));
end

nu=ceil(log2(n));%编码
code=zeros(length(a),nu);
for i=1:length(a)
    for j=nu:-1:0
        if (fix(b_quan(i)/(2^j))==1)
            code(i,(nu-j))=1;
            b_quan(i)=b_quan(i)-2^j;
        end
    end
end

%A率非线性的逆运算

a_quan = zeros(size(a));
for i = 1:length(c_quan);
    if( abs(c_quan(i)) <= 1/(1+log(A)) )
        a_quan(i) = (1+log(A))*c_quan(i)/A;
    end
    if(abs(c_quan(i))>1/A)
        a_quan(i) = sign(c_quan(i))*exp((1+log(A))*abs(c_quan(i))-1)/A;
    end
end
a_quan = a_quan*amax;
sqnr=20*log10(norm(a)/norm(a-a_quan));%求量化信噪比

% disp('量化信噪比')
% disp(sqnr)