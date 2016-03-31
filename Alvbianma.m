function [y,xx]=Alvbianma(x0)
x=x0';
Max=max(x);
x=x/Max;
z=sign(x);
x=abs(x);
for i=1:length(x)
    if(x(i)>=0&x(i)<1/64)
        y(i)=16*x(i);
        else if(x(i)>=1/64&x(i)<1/32)
            y(i)=8*x(i)+1/8;
        else if(x(i)>=1/32&x(i)<1/16)
            y(i)=4*x(i)+2/8;
        else if(x(i)>=1/16&x(i)<1/8)
            y(i)=2*x(i)+3/8;
        else if(x(i)>=1/8&x(i)<1/4)
            y(i)=x(i)+4/8;
        else if(x(i)>=1/4&x(i)<1/2)
            y(i)=(1/2)*x(i)+5/8;
        else if(x(i)>=1/2&x(i)<=1)
            y(i)=1/4*x(i)+6/8;
        end
      end
     end
    end
   end
  end
 end
end
    y=z.*y;
    
    f=zeros(length(y),8);
    z=sign(y);
    y=y*128;
    y=fix(y);
    y=abs(y);
    for i=1:length(y)
        if(y(i)==128)
            y(i)=127.9;
        end
    end
    for i=1:length(y)
        for j=6:-1:0
            f(i,8-j)=fix(y(i)/(2^j));
            y(i)=mod(y(i),(2^j));
        end
    end
    for i=1:length(y)
        if(z(i)==-1)
            f(i,1)=0;
        end
        if(z(i)==1)
            f(i,1)=1;
        end
    end
    
    
    
 for i=1:length(f)
     xx(1,i)=f(i,1)*(f(i,2)*64+f(i,3)*32+f(i,4)*16+f(i,5)*8+f(i,6)*4+f(i,7)*2+f(i,8));
 end
 xx=xx/128;
 xx=xx*Max;
            


