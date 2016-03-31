function [compressParams] = compressImage(Pic,L)
    
    imgSize = size(Pic);
    blocks = im2col(Pic,blokSizes,'distinct');                                     %divide the pic into blocks for KSVD
    blocks = double(blocks);
    CloumnMeans = mean(blocks);
    blocks = blocks - repmat(mean(blocks),prod(blokSizes),1);

    [PicName,PathName,~] = uigetfile('d:\m_sou_file\杨老师的遥感数据\*.*','Selcet a file');
    if(PicName == 0)
        return;
    end
    Dksvd = load([PathName,PicName]); 
    Dksvd = normc(Dksvd);
    
    params.data = blocks;
    params.memusage = 'high';
    params.initdict = Dksvd;
    
    
    params.Tdata = L;
    coeffice = SparseCode(params); 
    
    [row,cloumn,coe] = find(coeffice);
    coeSize = size(coeffice);

    coe =round(coe);
    index = find(coe~=0);
    coe = coe(index);

    Min = min(coe(:));
    Max = max(coe(:));
    coeliner = linerGary(coe,0,255)+1;

    [zippedCoe, infoCoe] = huffencode(coeliner);

    row = row(index);
    [zippedRow, infoRow] = huffencode(row);
    
    compressParams.zippedCoe = zippedCoe;
    compressParams.zippedRow = zippedRow;
    compressParams.infoCoe = infoCoe;
    compressParams.infoRow = infoRow;
    compressParams.CloumnMeans = CloumnMeans;
%     compressParams.cloumn = cloumn;
    compressParams.Min = Min;
    compressParams.Max = Max;
    compressParams.L = L;
    compressParams.Dksvd = Dksvd;
    compressParams.coeSize = coeSize;
    compressParams.imgSize = imgSize;
end
    