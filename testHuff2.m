
clc
clear
img =imread('lenna.bmp');
img = rgb2gray(img);
% img = img(ceil(end/4):ceil(end/4*2),ceil(end/4):ceil(end/4*2));
imgsize = size(img);
imshow(img);
img = double(img);
   
    
tic;
[r_encode,infor] = matrix2huff(img);
reImg = huff2matrix(r_encode,infor);

reImg = uint8(reImg);

figure(2);
imshow(reImg);
t(1) = toc;

tic;
[img_encoded, info] = huffencode(img(:));
img_decoed = huffdecode(img_encoded, info);
img_decoed = reshape(img_decoed,imgsize(1),imgsize(2));
img_decoed = uint8(img_decoed);

figure(3);
imshow(img_decoed);
t(2) = toc;

disp(num2str(t(1)));
disp(num2str(t(2)));

disp(sum(sum(img ~= reImg)));
disp(sum(sum(img ~= img_decoed)));

