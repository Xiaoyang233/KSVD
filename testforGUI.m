
Pic = imread('d:\m_sou_file\����ʦ��ң������\ȫɫͼ��\����\pan3.tif ');
ratio = 32;
[zippedCoe, infoCoe,zippedRow, infoRow,otherInfo] =  ksvdCompress(Pic,ratio);
[reconstKSVD] = ksvddecompress(zippedCoe, infoCoe,zippedRow, infoRow,otherInfo);

[PSNROutKsvd,~] = psnr(reconstKSVD,Pic);
imshow(reconstKSVD,[]);
title([' K-SVD (compression ratio: 32) ',num2str(PSNROutKsvd),'db ']);

