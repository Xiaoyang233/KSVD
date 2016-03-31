function [IOut,output,BppDct] = ImageDCT(Image,K,varargin)
%==========================================================================
%   P E R F O R M   D E N O I S I N G   U S I N G   O V E R C O M P L E T E 
%                        D C T    D I C T I O N A R Y
%==========================================================================
% function IOut = denoiseImageDCT(Image,sigma,bb,K)
% denoise an image by sparsely representing each block with the
% overcomplete DCT Dictionary, and averaging the represented parts.
% Detailed description can be found in "Image Denoising Via Sparse and Redundant
% representations over Learned Dictionaries", (appeared in the 
% IEEE Trans. on Image Processing, Vol. 15, no. 12, December 2006).
% ===================================================================
% INPUT ARGUMENTS : Image - the noisy image (gray-level scale)
%                   sigma - the s.d. of the noise (assume to be white Gaussian).
%                   K - the number of atoms in the representing dictionary.
%    Optional argumeters:              
%                  'blockSize' - the size of the blocks the algorithm
%                       works. All blocks are squares, therefore the given
%                       parameter should be one number (width or height).
%                       Default value: 8.
%                  'errorFactor' - a factor that multiplies sigma in order
%                       to set the allowed representation error. In the
%                       experiments presented in the paper, it was set to 1.15
%                       (which is also the default value here).
%                  'maxBlocksToConsider' - maximal number of blocks that
%                       can be processed. This number is dependent on the memory
%                       capabilities of the machine, and performances?
%                       considerations. If the number of available blocks in the
%                       image is larger than 'maxBlocksToConsider', the sliding
%                       distance between the blocks increases. The default value
%                       is: 250000.
%                  'slidingFactor' - the sliding distance between processed
%                       blocks. Default value is 1. However, if the image is
%                       large, this number increases automatically (because of
%                       memory requirements). Larger values result faster
%                       performances (because of fewer processed blocks).
%                  'waitBarOn' - can be set to either 1 or 0. If
%                       waitBarOn==1 a waitbar, presenting the progress of the
%                       algorithm will be displayed.
% OUTPUT ARGUMENTS : IOut - a 2-dimensional array in the same size of the
%                   input image, that contains the cleaned image.
%                    output - a struct that contains that following field:
%                       D - the dictionary used for denoising
% =========================================================================
%  初始化一个扁矩阵,其中，
%  K表示字典的大小,比如K取1024；
%  bb为字典原子维数开平方，如你要产生256为的，则设置bb=16；
%
%
%
% =========================================================================
Reduce_DC = 0;
[NN1,NN2] = size(Image);
C = 1.15;
waitBarOn = 1;
maxBlocksToConsider = 260000;
bb = 8;
for argI = 1:2:length(varargin)
    if (strcmp(varargin{argI}, 'slidingFactor'))
        slidingDis = varargin{argI+1};
    end
    if (strcmp(varargin{argI}, 'errorFactor'))
        C = varargin{argI+1};
    end
    if (strcmp(varargin{argI}, 'maxBlocksToConsider'))
        maxBlocksToConsider = varargin{argI+1};
    end
    if (strcmp(varargin{argI}, 'blockSize'))
        bb = varargin{argI+1};
    end
    if (strcmp(varargin{argI}, 'waitBarOn'))
        waitBarOn = varargin{argI+1};
    end
    if (strcmp(varargin{argI}, 'L'))
        L = varargin{argI+1};
    end
end
% Create an initial dictionary from the DCT frame
Pn=ceil(sqrt(K));
DCT=zeros(bb,Pn);
for k=0:1:Pn-1,
    V=cos([0:1:bb-1]'*k*pi/Pn);
    if k>0, V=V-mean(V); end;
    DCT(:,k+1)=V/norm(V);
end;
DCT=kron(DCT,DCT);
% while (prod(floor((size(Image)-bb)/slidingDis)+1)>maxBlocksToConsider)
%     slidingDis = slidingDis+1;
% end


% blocks = im2patch(Image, [bb,bb], [1,1]) ;
blocks = im2col(Image,[bb,bb],'distinct');

% % if (waitBarOn)
% %     counterForWaitBar = size(blocks,2);
% %     h = waitbar(0,'DCT compression In Process ...');
% % end
% go with jumps of 10000
% for jj = 1:10000:size(blocks,2)
% %     if (waitBarOn)
% %         waitbar(jj/counterForWaitBar);
% %     end
%     jumpSize = min(jj+10000-1,size(blocks,2)); % （jj+10000-1）为下一次循环所到达的列数
    if (Reduce_DC)
        vecOfMeans = mean(blocks(:,jj:jumpSize)); %每列求均值，得到行向量
        blocks(:,jj:jumpSize) = blocks(:,jj:jumpSize) - repmat(vecOfMeans,size(blocks,1),1);%｛blocks(:,jj:jumpSize)｝减去每列对应的均值
    end
% %     Coefs = OMPerr(DCT,blocks(:,jj:jumpSize),errT);
%      Coefs = OMP(DCT,blocks(:,jj:jumpSize),L);
      Coefs = OMP(DCT,blocks,L);
%=======================================================================
% encode part
[row,cloumn,coe] = find(Coefs);
coe1 =round(coe);
index = find(coe1~=0);
coe1 = coe1(index);

[totalNum,tempNum]=checkNum(coe1);  

[zippedCoe, infoCoe] = huffencode(coe1);
coe1 = huffdecode(zippedCoe, infoCoe);
coe1 = double(coe1);

row1 = row(index);
cloumn1 = cloumn(index);

[zippedRow, infoRow] = huffencode(row1);
row1 = huffdecode(zippedRow, infoRow);

[zippedClo, infoClo] = huffencode(cloumn1);
cloumn1 = huffdecode(zippedClo, infoClo);

[mm,nn]=size(Coefs);
Coefs1 = sparse(row1,cloumn1,coe1,mm,nn);

reBlock = DCT*Coefs1;

totalCodeNo=infoCoe.totalCodeNo+infoRow.totalCodeNo+infoClo.totalCodeNo;
bpp=totalCodeNo/(NN1*NN2);

% if (Reduce_DC)
%         blocks(:,jj:jumpSize)= DCT*Coefs + ones(size(blocks,1),1) * vecOfMeans;
%     else
% %         blocks(:,jj:jumpSize)= DCT*Coefs ;
%          blocks = DCT*Coefs ;
%     end
% end
% % if (waitBarOn)
% %     close(h);
% % end
BppDct=bpp;
% IOut = blocks;
IOut = reBlock;
output.D = DCT;