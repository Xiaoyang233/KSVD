

clc
clear
close
addpath('d:\m_sou_file\高光谱数据\压缩图像\');
prompt = {'the band compressed','the number of atoms perserved when compression'};
lineno = 1;
defans = {'9','8'};
nameTitle = 'input the parameter';
answer = inputdlg(prompt,nameTitle,lineno,defans);

for cellnum=1:length(answer)
        parameter(cellnum) = str2double(answer{cellnum});
end
band = parameter(1);
Ratio = parameter(2);

[file, pathname, filterindex] = uigetfile('d:\m_sou_file\高光谱数据\压缩图像\*.*', 'Pick a file');
if file == 0
    return;
end
[pathstr, fname, ext] = fileparts([pathname file]);
inputdir = [pathstr '\'];

% Pic = reshape(pixels(band,:),145,145);
Pic = imread([fname,ext]);
Pic_8bit = Pic;

% fileName = ['IndianaBand',num2str(band),'_',num2str(Ratio),'.jp2'];
fileName = ['jasper','_CR',num2str(Ratio),'_',fname,'.jp2'];
% Pic = uint16(Pic);
% Pic_8bit = uint8(linerGary(Pic),0,255);

jp2write(Pic,fileName,'rate',1/Ratio,'bitdepth',8);

figure(1);
imshow(Pic,[]);
title('the pervious image');

imwrite(Pic,['d:\m_sou_file\高光谱数据\压缩图像结果\',fname,'.bmp']);

JEPG2000Pic = jp2read(fileName);
% Pic_8bit = uint8(linerGary(Pic,0,255));
% JEPG2000Pic_8bit = uint8(linerGary(JEPG2000Pic,0,255)); 
% [PSNROutJEPG2000,~] = psnr(JEPG2000Pic_8bit,Pic_8bit);
[PSNROutJEPG2000,~] = psnr(JEPG2000Pic,Pic_8bit);

figure(2);
imshow(JEPG2000Pic,[]);
title(['the recovery image by JEPG2000 ',num2str(PSNROutJEPG2000),'db ']);

imwrite(JEPG2000Pic,['d:\m_sou_file\高光谱数据\压缩图像结果\',fileName(1:end-4),'.bmp']);

