function [reconstKSVD] = ksvddecompress(zippedCoe, infoCoe,zippedRow, infoRow,otherInfo)

    coeSize = otherInfo.coeSize;
    Min = otherInfo.Min;
    Max = otherInfo.Max;
    CloumnMeans = otherInfo.CloumnMeans;
    imgSize = otherInfo.imgSize;
    L = otherInfo.L;
    blokSize = 16;
    blokSizes=[blokSize,blokSize];
    load trainDicPan16_16;
    
    re_coe = huffdecode(zippedCoe, infoCoe);
    re_coe = double(re_coe);
    re_coe = linerGary(re_coe-1,Min,Max);
    re_coe = round(re_coe);

    re_row = huffdecode(zippedRow, infoRow);
    
    temp = 1:coeSize(2);
    temp = repmat(temp,L,1);
    cloumn = temp(:);

    coeffice = sparse(re_row,cloumn,re_coe,coeSize(1),coeSize(2));

    reBlock = Dksvd*coeffice;
    reBlock = reBlock + repmat(CloumnMeans,prod(blokSizes),1);

    if(isfield(otherInfo,'numblokClo'))
        numblokClo = otherInfo.numblokClo;
        j = 1;
        reconstKSVD = zeros(imgSize);
        for i = 1:imgSize(3)
            reconstKSVD(:,:,i) = col2im(reBlock(:,j:j + numblokClo-1),blokSizes,[imgSize(1),imgSize(2)],'distinct');
            j = j + numblokClo;
        end
    else
        reconstKSVD = col2im(reBlock,blokSizes,imgSize,'distinct');                  %combine blocks into pic
    end

    reconstKSVD = uint8(reconstKSVD);
end
            