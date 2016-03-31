function[zippedCoe, infoCoe,zippedRow, infoRow,otherInfo] =  ksvdCompress(Pic,ratio)
    load trainDicPan16_16;
    blokSize = 16;
    blokSizes=[blokSize,blokSize];

    imgSize = size(Pic);
    dim = length(imgSize);
    %% trun to block
    
    if(dim==3)
        j = 1;
        numblokClo = ceil(imgSize(1)/blokSize)*ceil(imgSize(2)/blokSize);
        blocks = zeros(prod(blokSizes),numblokClo*imgSize(3));
        for i = 1:imgSize(3)
            blocks(:,j:j + numblokClo - 1) = im2col(Pic(:,:,i),blokSizes,'distinct'); 
            j = j + numblokClo;
        end
    else
        blocks = im2col(Pic,blokSizes,'distinct');                                     %divide the pic into blocks for KSVD
    end
    
    blocks = double(blocks);
    CloumnMeans = mean(blocks);
    blocks = blocks - repmat(mean(blocks),prod(blokSizes),1);
    Idx = find(blocks == 0);
    blocks(Idx) = rand(length(Idx),1)*1e-12;
        zeroClo = checkSingal(blocks);

    %% compress with sparse code

    params.data = blocks;
    params.memusage = 'high';
    params.initdict = Dksvd;
    Dksvd = normc(Dksvd);
    L = prod(blokSizes)/ratio;
    params.Tdata = L;
                                                %sparse coding
    coeffice = SparseCode(params);  
    coeSize = size(coeffice);
    [row,cloumn,coe] = find(coeffice);
    [tem, num] = findNonzero(cloumn);
    len = findgap(tem);
    temp = 1:coeSize(2);
    for i = 1:length(zeroClo)
        tempIdx = find(temp == zeroClo(i));
%         if(tempIdx(end)<length(tempIdx))
            temp = [temp(1:tempIdx(1)-1),temp(tempIdx(end)+1:end)];
%         else
%             temp = [temp(1:tempIdx(1)-1),temp(tempIdx(end)+1)];
    end
    temp = repmat(temp,L,1);
    cloumn1 = temp(:);

    coe =round(coe);
    index = find(coe~=0);
    coe = coe(index);

    Min = min(coe(:));
    Max = max(coe(:));
    coe21 = linerGary(coe,0,255)+1;

    [zippedCoe, infoCoe] = huffencode(coe21);

    row = row(index);
    [zippedRow, infoRow] = huffencode(row);
    %% other infomation for decode
    otherInfo.Min = Min;
    otherInfo.Max = Max;
    otherInfo.CloumnMeans = CloumnMeans;
    otherInfo.imgSize = imgSize;
    if(exist('numblokClo','var'))
         otherInfo.numblokClo = numblokClo;
    end
    otherInfo.coeSize = coeSize;
    otherInfo.L = L;
end
function[tem, num] = findNonzero(vector)
    tem = unique(vector);
    num = zeros(length(tem),1);
    for i = 1:length(tem)
        num(i) = length(find(vector==tem(i)));
    end
end
function[len] = findgap(vector)
Max = max(vector(:));
len = zeros(length(vector),1);
    for i = 1:Max
        len(i) = length(find(vector==i));
    end
end
function [zeroClo] = checkSingal(singal)
    temp = singal>1e-12;
    temp = sum(temp);
    zeroClo = find(temp==0);
end