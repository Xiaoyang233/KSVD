clear all;
close all;
clc;
% i=uint8(0);
bppKSVD=[3.2238 1.6086 0.73181];
 psnrKSVD=[35.6444 31.2354 26.997];
 bppDCT=[3.0729,1.509,0.69739];
 psnrDCT=[33.2268,28.3095,24.6166];
 string=['   8±Ά' 
     '  16±Ά' 
     '  32±Ά'];
 figure(1);
 axis([0.5 3.6 25 36]);
 set(gca,'xtick',0.5:0.1:3.5);
 set(gca,'ytick',25:0.5:36);
  axis equal;
 plot(bppKSVD,psnrKSVD,'r-s','linewidth',2,'MarkerEdgeColor','y','MarkerFaceColor','b');
 hold on;
 plot(bppDCT,psnrDCT,'b-*','linewidth',2,'MarkerEdgeColor','r','MarkerFaceColor','b');
% plot(bppKSVD,psnrKSVD,bppDCT,psnrDCT,'r-s','r-s','linewidth',2,'MarkerEdgeColor','y','MarkerFaceColor','b');
 xlabel('bpp','FontSize',12, 'FontName','Times New Roman','fontweight','bold');
 ylabel('psnr','FontSize',12, 'FontName','Times New Roman','fontweight','bold');
 legend('KSVD','DCT','Upper bound','Analytical',2);
 legend('boxoff');
set(legend,'fontname','±κΏ¬Με');
set(legend,'fontweight','bold');
hold off;
 for i=1:3   
    text(bppKSVD(i),psnrKSVD(i),string(i,:));
    text(bppDCT(i),psnrDCT(i),string(i,:));
 end