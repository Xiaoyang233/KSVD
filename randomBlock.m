function[selectBlock] = randomBlock(Img,blkSize,selectNum)
% randomly select a squre in the image
% note: the size of squre must be even
% input:
% Img: image
% blkSize: 2 dimension vector decripting the size of squre
% selectNum: the num of squre
% output:
% selectBlock: the squre in the image

	ImgSize = size(Img);
	distance = ceil(blkSize/2);
	idxMax = zeros(ImgSize(1),ImgSize(2));
	idxMax(distance(1):ImgSize(1)-distance(1)+1,...
        distance(2):ImgSize(2)-distance(2)+1) = 1;
	tempIdx = find(idxMax(:));
	randNum = randperm(length(tempIdx),selectNum);
	randIdx = tempIdx(randNum);
	[row,col] = ind2sub(size(idxMax),randIdx);
	selectBlock = cell(1,selectNum);
    rowIdx = cell(1,selectNum);
    colIdx = cell(1,selectNum);
	for i = 1:selectNum
		rowIdx{i} = row(i) - distance(1):row(i) + distance(1)-1;
		colIdx{i} = col(i) - distance(2):col(i) + distance(2)-1;

		if(row(i) == distance(1))
			rowIdx{i} = rowIdx{i} + 1;
		end
		if(col(i) == distance(2))
			colIdx{i} = colIdx{i} + 1;
		end
	selectBlock{i} = Img(rowIdx{i},colIdx{i},:);
	end
end