function testHuff
    clc
    clear
    img =imread('lenna.bmp');
    img = rgb2gray(img);
    img = img(ceil(end/4):ceil(end/4*2),ceil(end/4):ceil(end/4*2));
    imgsize = size(img);
    imshow(img);
    
 	tic;
 	[template,Pro] = frequency(img(:));
 	dict = huffmandict(template, Pro);

 	img_encoded = huffmanenco(img(:),dict);
 	img_decoed = huffmandeco(img_encoded, dict);
    
    img_decoed = reshape(img_decoed,imgsize(1),imgsize(2));
 	figure(2);
 	imshow(img_decoed);
    t(1) = toc;
    disp(sum(sum(img ~= img_decoed)));    
    
    tic;
    [img_encoded, info] = huffencode(img(:));
    img_decoed = huffdecode(img_encoded, info);
    img_decoed = reshape(img_decoed,imgsize(1),imgsize(2));
    img_decoed = uint8(img_decoed);
    
    figure(3);
 	imshow(img_decoed);
    t(2) = toc;
    disp(num2str(t(1)));
    disp(' ');
    disp(num2str(t(2)));
    disp(sum(sum(img ~= img_decoed)));
end