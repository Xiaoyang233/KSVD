
function [r_encode,infor] = matrix2huff(matrix)
	[row,col] = size(matrix);
	eob = max(matrix(:)) + 1;              		 % Create end-of-block symbol
	r = zeros(numel(matrix) + size(matrix, 2), 1);	% the worst condition is that 
													% all the nonzero coeffice last of which is at end of cloumn plus symbol 'eob'
	count= 0;
	for j = 1:col                       		% Process 1 block (col) at a time
	    i = max(find(matrix(:, j)));        	 % Find last non-zero element
	    if isempty(i)                   	% No nonzero block values
	       i = 0;
	    end
	    p = count + 1;
	    q = p + i;
	    r(p:q) = [matrix(1:i, j); eob];      	% Truncate trailing 0's, add EOB,
	    count = count + i + 1;         		 % and add to output vector
	end

	r((count + 1):end) = [];         		  % Delete unusued portion of r
	[r_encode,infor] = huffencode(r);
	infor.M_row = row;
	infor.M_col = col;
end