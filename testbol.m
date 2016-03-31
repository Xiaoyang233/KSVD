clear all;
close all;
clc;
load Dic;
NoPic = 1000;
blokSize=[8,8];
L=8;
K=512;
sigma = 25;
sliding = 1;
selSample = randperm(NoPic-3,1)+3;
path = uigetdir('D:\','select the path of test pic');
PicVector = dir(path);
PicStr=strcat(path,'\',PicVector(selSample).name)
testPicRGB = imread(PicStr);                                                            %read the test pic
testPic = rgb2gray(testPicRGB)                                                      %RGB mode to Gray mode 
testPic = double(testPic);
% blocks = im2col(testPic,blokSize,'distinct');                                     %divide the pic into blocks
% blocks = im2col(testPic,blokSize,'sliding');
% blocks =double(blocks);
% reblocks = col2im(blocks,blokSize,size(testPic),'sliding');
blocks = im2patch(testPic, blokSize, [1,1]);
reblocks = patch2im(blocks, size(testPic), blokSize, [1,1]);
imshow(reblocks,[]);