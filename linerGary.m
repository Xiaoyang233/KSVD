function Y = linerGary(X,Ymin,Ymax)
% min-ͼ��Ҷ���Сֵ
% max-ͼ��Ҷ����ֵ
% X-����ͼ��
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