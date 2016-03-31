function Y = linerGary(X,Ymin,Ymax)
% min-图像灰度最小值
% max-图像灰度最大值
% X-输入图像
    X = double(X);
    Xmax = max(X(:));
    Xmin = min(X(:));
    if(Xmax==Xmin)
        error(message('the elements in the image cannot be all equal!'));
    end
    Y = (X(:)-Xmin)/(Xmax-Xmin)*(Ymax-Ymin)+Ymin;
    if length(size(X))==2
        [r,c] = size(X);
        Y = reshape(Y,r,c);
    else
        [r,c,d] = size(X);
        Y = reshape(Y,r,c,d);
    end
end