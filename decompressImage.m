function [reconstKSVD] = decompressImage(compressParams)

    zippedCoe = compressParams.zippedCoe;
    zippedRow = compressParams.zippedRow;
    infoCoe = compressParams.infoCoe;
    infoRow = compressParams.infoRow;
    CloumnMeans = compressParams.CloumnMeans;
%     cloumn = compressParams.cloumn;
    Min = compressParams.Min;
    Max = compressParams.Max;
    L = compressParams.L;
    Dksvd = compressParams.Dksvd;
    coeSize = compressParams.coeSize;
    imgSize = compressParams.imgSize;
    
    coe = huffdecode(zippedCoe, infoCoe);
    coe = double(coe);
    coe = linerGary(coe-1,Min,Max);
    coe = round(coe);
    
    row = huffdecode(zippedRow, infoRow);
    
    temp = 1:coeSize(2);
    temp = repmat(temp,L,1);
    cloumn = temp(:);
    
    coeffice = sparse(row,cloumn,coe,coeSize(1),coeSize(2));
    
    reBlock = Dksvd*coeffice;
    reBlock = reBlock + repmat(CloumnMeans,prod(blokSizes),1); 
    
    reconstKSVD = col2im(reBlock,blokSizes,imgSize,'distinct');                  %combine blocks into pic
    reconstKSVD = uint8(reconstKSVD);
    
end