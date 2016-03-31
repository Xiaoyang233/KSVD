function [matrix] = huff2matrix(r_encode,infor)

	r_decode = huffdecode(r_encode,infor);
    eob = max(r_decode(:));
	matrix = zeros(infor.M_row,infor.M_col);   k = 1;           % Form block columns by copying

	for j = 1:infor.M_col                          % successive values from x into
	   for i = 1:infor.M_row + 1                      % columns of matrix, while changing
	      if r_decode(k) == eob                  % to the next column whenever
	         k = k + 1;   break;          % an EOB symbol is found.
	      else
	         matrix(i, j) = r_decode(k);
	         k = k + 1;
	      end
	   end
	end
	
end