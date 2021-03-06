function [output] = ImageKSVD(blkMatrix,K,varargin)
%==========================================================================
%   P E R F O R M   D E N O I S I N G   U S I N G   A  D I C T  I O N A R Y
%                  T R A I N E D   O N      I M A G E  BLOCKSS
%==========================================================================
% function IOut = denoiseImageKSVD(Image,sigma,K,varargin)
% denoise an image by sparsely representing each block with the
% already overcomplete trained Dictionary, and averaging the represented parts.
% Detailed description can be found in "Image Denoising Via Sparse and Redundant
% representations over Learned Dictionaries", (appeared in the 
% IEEE Trans. on Image Processing, Vol. 15, no. 12, December 2006).
% This function may take some time to process. Possible factor that effect
% the processing time are:
%  1. number of KSVD iterations - the default number of iterations is 10.
%  However, fewer iterations may, in most cases, result an acceleration in
%  the process, without effecting  the result too much. Therefore, when
%  required, this parameter may be re-set.
%  2. maxBlocksToConsider - The maximal number of blocks to train on. If this 
%  number is larger the number of blocks in the image, random blocks
%  from the image will be selected for training. 
% ===================================================================
% INPUT ARGUMENTS : Image - the noisy image (gray-level scale)
%                   sigma - the s.d. of the noise (assume to be white Gaussian).
%                   K - the number of atoms in the trained dictionary.
%    Optional arguments:              
%                  'blockSize' - the size of the blocks the algorithm
%                       works. All blocks are squares, therefore the given
%                       parameter should be one number (width or height).
%                       Default value: 8.
%                       'errorFactor' - a factor that multiplies sigma in order
%                       to set the allowed representation error. In the
%                       experiments presented in the paper, it was set to 1.15
%                       (which is also the default  value here).
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
%                  'numKSVDIters' - the number of KSVD iterations processed
%                       blocks from the noisy image. If the number of
%                       blocks in the image is larger than this number,
%                       random blocks from all available blocks will be
%                       selected. The default value for this parameter is:
%                       10 if sigma > 5, and 5 otherwise.
%                  'maxNumBlocksToTrainOn' - the maximal number of blocks
%                       to train on. The default value for this parameter is
%                       65000. However, it might not be enough for very large
%                       images
%                  'displayFlag' - if this flag is switched on,
%                       announcement after finishing each iteration will appear,
%                       as also a measure concerning the progress of the
%                       algorithm (the average number of required coefficients
%                       for representation). The default value is 1 (on).
%                  'waitBarOn' - can be set to either 1 or 0. If
%                       waitBarOn==1 a waitbar, presenting the progress of the
%                       algorithm will be displayed.
% OUTPUT ARGUMENTS : Iout - a 2-dimensional array in the same size of the
%                       input image, that contains the cleaned image.
%                    output.D - the trained dictionary.
% =========================================================================

reduceDC = 0;
% [NN1,NN2] = size(blkMatrix);
waitBarOn = 1;
param.numIteration = 25; %ksvd��������

bb = 8;

displayFlag = 1;
% maxNumBlocksToTrainOn = 65000;
% maxBlocksToConsider = 260000;
% slidingDis = 1;

param.K = K;
param.errorFlag = 0;
param.L=8;
param.preserveDCAtom = 1;

for argI = 1:2:length(varargin)
%     if (strcmp(varargin{argI}, 'slidingFactor'))
%         slidingDis = varargin{argI+1};
%     end
%     if (strcmp(varargin{argI}, 'maxBlocksToConsider'))
%         maxBlocksToConsider = varargin{argI+1};
%     end
    if (strcmp(varargin{argI}, 'numKSVDIters'))
        param.numIteration = varargin{argI+1};
    end
    if (strcmp(varargin{argI}, 'blockSize'))
        bb = varargin{argI+1};
    end
%     if (strcmp(varargin{argI}, 'maxNumBlocksToTrainOn'))
%         maxNumBlocksToTrainOn = varargin{argI+1};
%     end
    if (strcmp(varargin{argI}, 'displayFlag'))
        displayFlag = varargin{argI+1};
    end
    if (strcmp(varargin{argI}, 'waitBarOn'))
        waitBarOn = varargin{argI+1};
    end
    if (strcmp(varargin{argI},'reduceDC'))
        reduceDC = varargin{argI+1};
    end
    if (strcmp(varargin{argI}, 'ModeOMP'))
        param.errorFlag = varargin{argI+1};
    end
    if (strcmp(varargin{argI}, 'FirstAtomPre'))
        param.preserveDCAtom = varargin{argI+1};
    end
    if (strcmp(varargin{argI}, 'NO_atoms'))
        param.L = varargin{argI+1};
    end
    if (strcmp(varargin{argI}, 'errorGoal'))
        param.errorGoal = varargin{argI+1};
    end
end

% first, train a dictionary on blocks from the noisy image

% maxNumBlocksToTrainOn:  at this time,the prameter 
%choose 'sliding' and the slidingFactor is '1',thus prod([NN1,NN2]-bb+1
%is maxNumBlocksToTrainOn

%if prod([NN1,NN2]-bb+1)> maxNumBlocksToTrainOn��slidingFactor=1

% if(prod([NN1,NN2]-bb+1)> maxNumBlocksToTrainOn)
%     randPermutation =  randperm(prod([NN1,NN2]-bb+1));
%     selectedBlocks = randPermutation(1:maxNumBlocksToTrainOn);
% 
%     blkMatrix = zeros(bb^2,maxNumBlocksToTrainOn);
%     for i = 1:maxNumBlocksToTrainOn
%         [row,col] = ind2sub(size(Image)-bb+1,selectedBlocks(i));
%         currBlock = Image(row:row+bb-1,col:col+bb-1);
%         blkMatrix(:,i) = currBlock(:);
%     end
% else
%     blkMatrix = im2col(Image,[bb,bb],'sliding');
% end

Pn=ceil(sqrt(K));
DCT=zeros(bb,Pn);
for k=0:1:Pn-1,
    V=cos([0:1:bb-1]'*k*pi/Pn);
    if k>0, V=V-mean(V); end;
    DCT(:,k+1)=V/norm(V);
end;
DCT=kron(DCT,DCT);

param.initialDictionary = DCT(:,1:param.K );
param.InitializationMethod =  'GivenMatrix';

if (reduceDC)
    vecOfMeans = mean(blkMatrix);
    blkMatrix = blkMatrix-ones(size(blkMatrix,1),1)*vecOfMeans;
end

if (waitBarOn)
    counterForWaitBar = param.numIteration+1;
    h = waitbar(0,'Traning In Process ...');
    param.waitBarHandle = h;
    param.counterForWaitBar = counterForWaitBar;
end


param.displayProgress = displayFlag;
[Dictionary,output] = KSVD(blkMatrix,param);
output.D = Dictionary;

if (displayFlag)
    disp('finished Trainning dictionary');
end
if (waitBarOn)
    close(h);
end