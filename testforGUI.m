
Pic = imread('d:\m_sou_file\杨老师的遥感数据\全色图像\西安\pan3.tif ');
ratio = 32;
[zippedCoe, infoCoe,zippedRow, infoRow,otherInfo] =  ksvdCompress(Pic,ratio);
[reconstKSVD] = ksvddecompress(zippedCoe, infoCoe,zippedRow, infoRow,otherInfo);

[PSNROutKsvd,~] = psnr(reconstKSVD,Pic);
imshow(reconstKSVD,[]);
title([' K-SVD (compression ratio: 32) ',num2str(PSNROutKsvd),'db ']);

