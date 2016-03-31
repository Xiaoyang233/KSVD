function [template,Pro] = frequency(data)

    % data: vector
    % Pro: the probility of each one in data
	
    if(isvector(data))
		template = unique(data);
    	UniqueNum = length(template);
        totalNum = length(data);
    	Pro = zeros(1,UniqueNum);
        
            for i = 1:UniqueNum
                num = length(find(data == template(i)));
                Pro(i) = num/totalNum;
            end
    end
end