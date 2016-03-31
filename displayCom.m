function [] = displayCom(Pic,reconstKSVD)

    [PSNROutKsvd,~] = psnr(reconstKSVD,Pic);

    figure(1);
    subplot(1,3,1);
    imshow(Pic,[]);
    title('the pervious image');

    subplot(1,3,2);
    JEPG2000Pic = jp2read(fileName);
    imshow(JEPG2000Pic,[]);

    [PSNROutJEPG2000,~] = psnr(JEPG2000Pic,Pic);
    title([' JEPG2000 (compression ratio: 8) ',num2str(PSNROutJEPG2000),'db ']);

    subplot(1,3,3);

    imshow(reconstKSVD,[]);
    % title(['the recovery image by K-SVD ',num2str(PSNROutKsvd),'db ',num2str(bpp),'bpp']);
    title([' K-SVD (compression ratio: 8) ',num2str(PSNROutKsvd),'db ']);
end