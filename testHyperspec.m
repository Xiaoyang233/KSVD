clc;
clear all;
close all;
load Salinas-A;
load Botswana_data;
load IndianaPines;
figure(1);
% Indiana = multibandread('IndianaPines.mat', [145,145,200], 'uint8=>uint8',...
%                      0, 'bil', 'ieee-le', {'Band','Direct',1});
% imshow(Indiana);
% xlabel('Image Indiana Space Imaging');
% 
% figure(2);
% Botswana = multibandread('Botswana_data.mat', [145,145,200], 'uint8=>uint8',...
%                      0, 'bil', 'ieee-le', {'Band','Range',[1 3]});
% imshow(Botswana);
% xlabel('Image Botswana Space Imaging');
% 
% figure(3);
% Salinas = multibandread('Salinas-A.mat', [145,145,1], 'uint8=>uint8',...
%                      0, 'bil', 'ieee-le', {'Band','Range',[1 3]});
% imshow(Salinas);
% xlabel('Image Salinas Space Imaging');

Max=max(max(pixels));
Min=min(min(pixels));
R=124;
G=17;
B=35;
tempR =reshape(pixels(R,:),[145,145]);
tempG =reshape(pixels(G,:),[145,145]);
tempB =reshape(pixels(B,:),[145,145]);
tempR =normlization(tempR,Max,Min);
tempG =normlization(tempG,Max,Min);
tempB =normlization(tempB,Max,Min);

temp=cat(3,tempR,tempG,tempB);
imshow(temp,[]);