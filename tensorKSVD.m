% Find a s-sparse representation of a tensor Y, i.e. 
% Y = X \times_1 D_1 \times_2 D_2...\times_N D_N
% where X(i_1,i_2,...,i_N)=0 for all i_n \notin I_n and |I_n|=s(n)


function [X,D] = tensorKSVD(D,Y,s,epsilon,numIteration)
    [I] = size(Y);
    N = size(I,2);  % the dimension of Y

    for n = 1:N
        M(n) = size(D{n},2);     % D is cell，each cell stores a dictionary;
                                 % each element in M store the number of atom in the coresponding dictionary
    end

    norma = norm(reshape(Y,[I(1),prod(I)/I(1)]),'fro');     % the sum of squres of elements in Y, and then squre root it

    % Main loop where selected indices in each mode are found
    %Dred = D;
    reconstError = zeros(1,numIteration);
    for Iteration = 1:numIteration
        
        R = Y; % initial residual
        Ind = cell(1,N); % set of indices to select for each mode
        L = cell(1,N); % Triangular matrices used in Cholesky factorization of each mode
        Dsub = cell(1,N); % subset of selected atoms per each mode

        % auxiliar variables for searching the maximum of a tensor
        v = cell(1,N); 
        add = cell(1,N); 

        ni = zeros(1,N); % number of selected indices in each mode

        % Initialization
        for n = 1:N
            Ind{n} = [];
            L{n} = 1;
            ni(n) = 0;
        end
        X = zeros(M); % coefficients
        cond = 0; % this condition is true when ni(n)=s(n) for all n (all nonzero coefficients where computed)
        posi = zeros(1,N); % additional index to add


        % multiway cross correlation
        B = double(ttensor(tensor(Y),transp(D)));       % resiual is initialy set as Y，
                                                        % B is multiway cross correlation between dictionary and resiual

        condchange = 1;
        error = Inf;
        
        allInd = cell(1,N);
    
        while ((condchange && (~cond) && (error > epsilon)))    
            niant = ni;
            proj = abs(double(ttensor(tensor(R),transp(D)))); % multiway correlation between dictionary and residual
                                                              % the diffence between 'B' and 'proj' is that 'proj' find the absolute value
            % Proj{1} = proj;
            % [Proj{2},t{1}] = max(abs(Proj{1}));
            % [Proj{3},t{2}] = max(abs(Proj{2}));
            % [Proj{4},t{3}] = max(abs(Proj{3}));
        %     三维情况：寻找最?索引，v{1}:张量中每列所对应的最?索引，把每列的最?组成新的矩阵
        %     v{2}:新的矩阵对应每行的最?索引，每行的???值组成新的矢?
        %     v{3}:新的矢量对应的最?索引
            for n = 1:N
                [proj,v{n}] = max(abs(proj));
                add{n} = 1;
            end
        %     确定整个张量中最?的索?v{3}???值第三维索引，v{2}[1,1,v{3}]???值第二维索引?
        %     v{2}是矢量，v{3}个元素存储的是求??前它???行索引；
        %     v{1}[1,v{2},v{3}]???值的第一维索?
            posi(N) = v{N};
            add{N} = posi(N);
            for n = N-1:-1:1        
                posi(n) = v{n}(add{:});
                add{n} = posi(n);
            end

        %     proj = reshape(proj,[prod(M),1]);
        %      
        %     [maxVal,pos]=max(proj);
        %     pos=pos(1);
        %     
        %     posi = multi_index(M,pos);
        %     for n = 1:N
        %         add{n} = posi(n);
        %     end            
            
            for n = 1:N
                allInd{n} = [allInd{n},posi(n)];

                if ((~ismember(posi(n),Ind{n})) && (ni(n) < s(n)))
                Ind{n} = [Ind{n},posi(n)];
                Dsub{n} = [Dsub{n},D{n}(:,posi(n))];
                ni(n) = ni(n) + 1;
                    if ni(n) > 1
                        w = L{n}\((D{n}(:,Ind{n}(1:ni(n) - 1)))'*D{n}(:,Ind{n}(ni(n))));
                        L{n} = [L{n}, zeros(ni(n)-1,1); w', sqrt(1 - w'*w)];
                    end         
                end
            end
            
            Z = B(Ind{:});
            
            %% Aca esta el problema, comparar con el codigo
            for n = 1: N
                Z = L{n}\reshape(permute(Z,permorder(n,N)),[ni(n),prod(ni)/ni(n)]);
                Z = L{n}'\Z;
                Z = permute(reshape(Z,ni(permorder(n,N))),permorderinv(n,N));
            end
             Z = tensor(Z,ni);                      % 求出的系?
            

        %   
        %      Z = L{1}\reshape(permute(B(Ind{1},Ind{2},Ind{3}),permorder(1,3)),[ni(1),ni(2)*ni(3)]);
        %      X1 = L{1}'\Z;
        %      X1 = permute(reshape(X1,[ni(1),ni(2),ni(3)]),permorderinv(1,3));
        %      
        %      Z = L{2}\reshape(permute(X1,permorder(2,3)),[ni(2),ni(1)*ni(3)]);
        %      X2 = L{2}'\Z;
        %      X2 = permute(reshape(X2,[ni(2),ni(1),ni(3)]),permorderinv(2,3));   
        %      
        %      Z = L{3}\reshape(permute(X2,permorder(3,3)),[ni(3),ni(1)*ni(2)]);     
        %      X3 = L{3}'\Z;
        %      Z = permute(reshape(X3,[ni(3),ni(1),ni(2)]),permorderinv(3,3));
        %      
        %      Z = tensor(Z,[ni(1),ni(2),ni(3)]);
            
            
            
            %Xsub1 = double(ttensor(tensor(Y),{pinv(D1(:,I1)),pinv(D2(:,I2)),pinv(D3(:,I3))}));
            
            %norm(reshape(Xsub-Xsub1,[ni1,ni2*ni3]),'fro')
           
            
            R = Y - double(ttensor(Z,Dsub));                            % 残差 = 原始信号 - 块稀疏系?对应的子字典
            error = norm(reshape(R,[I(1),prod(I)/I(1)]),'fro')/norma;
        %     disp([num2str(ni),'  ', num2str(error)])          %为了节省时间，不显示
            
            cond = 1;
            for n =1:N
                cond = cond && (ni(n) == s(n));
            end
            
            condchange = sum(ni-niant);
            
        %     for n = 1:N
        %         Dred{n}(:,Ind{n}(ni(n)))=0;
        %     end
        end
    

        X(Ind{:}) = Z;
        X=sptensor(X);

        coeffNum = length(allInd{1});
        rPerm = randperm(coeffNum);
        Posi = zeros(1,N);
        for i = rPerm
            for j = 1:N
                Posi(j) = allInd{j}(i);
                % tempD{j} = D{j}(i);
            end

            X(Posi(1),Posi(2),Posi(3)) = 0;             % remain to change
            Error = Y - double(ttensor(tensor(X),D));
            Error = sptensor(Error);
            P = parafac_als(Error,1);

            for j = 1:N
                D{j}(:,Posi(j)) = P.U{j};
            end
            X(Posi(1),Posi(2),Posi(3)) = P.lambda;      % remain to change
        end

        reconstError(Iteration) = (norm(tensor(Y - double(full(ttensor(X,D)))))^2)/prod(size(Y));
        disp(['iteration ',num2str(Iteration),' :Error is ',num2str(reconstError(Iteration))]);
        plot(reconstError);
    end
end
function [D] = transp(D)
N = size(D,2);
for n = 1:N
    D{n} = D{n}';
end
end

function [v] = permorder(n,N)
% 把v中的前n个数循环向右移动1位，其他的N-n个数数保持不变；
v = 1:N;
if n == 1 
    return;
end
v(1) = n;
for m = 1:n-1
    v(m+1) = m; 
end
end

function [v] = permorderinv(n,N)
v = 1:N;
if n == 1 
    return;
end
v(1:n-1) = 2:n;
v(n) = 1;

end

function [indx] = multi_index(M,pos)
N = size(M,2);
indx = zeros(1,N);
for n = N:-1:2
    indx(n) = floor(((pos - 1)/prod(M(1:n-1)))) + 1;
    pos = pos - (indx(n)-1)*prod(M(1:n-1));
end
indx(1) = pos;

end
% function [betterDictionaryElement,X,NewVectorAdded] = I_findBetterDictionaryElement( Y,D,j,X,numCoefUsed)

% % X: core tensor
% % Y: signal tensor
% % D: dictionary
% if (length(who('numCoefUsed'))==0)
%     numCoefUsed = 1;
% end
% relevantDataIndices = find(X(j,:)); % the Y indices that uses the j'th D element.
% if (length(relevantDataIndices)<1) %(length(relevantDataIndices)==0)
%     ErrorMat = Y-D*X;
%     ErrorNormVec = sum(ErrorMat.^2);
%     [d,i] = max(ErrorNormVec);
%     betterDictionaryElement = Y(:,i);%ErrorMat(:,i); %
%     betterDictionaryElement = betterDictionaryElement./sqrt(betterDictionaryElement'*betterDictionaryElement);
%     betterDictionaryElement = betterDictionaryElement.*sign(betterDictionaryElement(1));
%     X(j,:) = 0;
%     NewVectorAdded = 1;
%     return;
% end

% NewVectorAdded = 0;
% tmpCoefMatrix = X(:,relevantDataIndices); 
% tmpCoefMatrix(j,:) = 0;% the coeffitients of the element we now improve are not relevant.
% errors = Y(:,relevantDataIndices) - D*tmpCoefMatrix); % vector of errors that we want to minimize with the new element
% % % the better dictionary element and the values of beta are found using svd.
% % % This is because we would like to minimize || errors - beta*element ||_F^2. 
% % % that is, to approximate the matrix 'errors' with a one-rank matrix. This
% % % is done using the largest singular value.
% [betterDictionaryElement,singularValue,betaVector] = svds(errors,1);
% X(j,relevantDataIndices) = singularValue*betaVector';% *signOfFirstElem
