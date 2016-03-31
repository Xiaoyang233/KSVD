% Find a s-sparse representation of a tensor Y, i.e. 
% Y = X \times_1 D_1 \times_2 D_2...\times_N D_N
% where X(i_1,i_2,...,i_N)=0 for all i_n \notin I_n and |I_n|=s(n)


function [X,D] = tensorKSVD1(D,Y,s,epsilon,numIteration)
    [I] = size(Y);
    N = size(I,2);  % the dimension of Y

    for n = 1:N
        M(n) = size(D{n},2);     % D is cellï¼Œeach cell stores a dictionary;
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
        B = double(ttensor(tensor(Y),transp(D)));       % resiual is initialy set as Yï¼?
                                                        % B is multiway cross correlation between dictionary and resiual

        condchange = 1;
        error = Inf;
        
        allInd = cell(1,N);
    
        while ((condchange && (~cond) && (error > epsilon)))    
            niant = ni;
            proj = abs(double(ttensor(tensor(R),transp(D)))); % multiway correlation between dictionary and residual
                                                              % the diffence between 'B' and 'proj' is that 'proj' find the absolute value

            for n = 1:N
                [proj,v{n}] = max(abs(proj));
                add{n} = 1;
            end

            posi(N) = v{N};
            add{N} = posi(N);
            for n = N-1:-1:1        
                posi(n) = v{n}(add{:});
                add{n} = posi(n);
            end
           
            
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
             Z = tensor(Z,ni);                      % æ±‚å�?ºçš„ç³�?
            

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
           
            
            R = Y - double(ttensor(Z,Dsub));                     % resiual = signal - the subdictionary for block sparse coeffice
            error = norm(reshape(R,[I(1),prod(I)/I(1)]),'fro')/norma;
            %     disp([num2str(ni),'  ', num2str(error)])          % no display for saving time
            
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
        
        Xmode_1 = reshape(X,M(1),M(2)*M(3));    
        for i = 1:M(1)
            if(ismember(i,Ind{1}))

                tempXmode_1 = Xmode_1;
                tempXmode_1(i,:) = 0;

                tempResult = D{1}*tempXmode_1;
                tempCore = reshape(tempResult,I(1),M(2),M(3));
                tempMode_2 = reshape(shiftdim(tempCore,1),M(2),I(1)*M(3));
                tempResult = D{2}*tempMode_2;
                tempCore = shiftdim(reshape(tempResult,I(2),M(3),I(1)),1);  % we have shifted 90 degree aforesaid, 
                                                                            % now shift 90 degree again , total 180 degree
                tempMode_3 = reshape(tempCore,M(3),I(1)*I(2));
                tempResult = D{3}*tempMode_3;
                tempCore = shiftdim(reshape(tempResult,I(3),I(1),I(2)),1);

                Error = Y - tensor(tempCore);
                Error = sptensor(Error);
                P = parafac_als(Error,1);
                D{1}(:,i) = P.U{1};
            end
        end

        clear Xmode_1 tempXmode_1 tempMode_2 tempMode_3 Error tempResult tempCore;

        Xmode_2 = reshape(shiftdim(X,1),M(2),M(1)*M(3));
        for i = 1:M(2)
            if(ismember(i,Ind{2}))

                tempXmode_2 = Xmode_2;
                tempXmode_2(i,:) = 0;

                tempResult = D{2}*tempXmode_2;
                tempCore = reshape(tempResult,I(2),M(3),M(1));

                tempMode_3 = reshape(shiftdim(tempCore,1),M(3),M(1)*I(2));
                tempResult = D{3}*tempMode_3;
                tempCore = reshape(tempResult,I(3),M(1),I(2));

                tempMode_1 = reshape(shiftdim(tempCore,1),M(1),I(2)*I(3));
                tempResult = D{1}*tempMode_1;
                tempCore = reshape(tempResult,I(1),I(2),I(3)); 

                Error = Y - tensor(tempCore);
                Error = sptensor(Error);
                P = parafac_als(Error,1);
                D{2}(:,i) = P.U{2};
            end
        end

        clear Xmode_2 tempXmode_2 tempMode_1 tempMode_3 Error tempResult tempCore;

        Xmode_3 = reshape(shiftdim(X,2),M(3),M(1)*M(2));
        for i = 1:M(3)
            if(ismember(i,Ind{3}))

                tempXmode_3 = Xmode_3;
                tempXmode_3(i,:) = 0;

                tempResult = D{3}*tempXmode_3;
                tempCore = reshape(tempResult,I(3),M(1),M(2));

                tempMode_1 = reshape(shiftdim(tempCore,1),M(1),M(2)*I(3));
                tempResult = D{1}*tempMode_1;
                tempCore = reshape(tempResult,I(1),M(2),I(3));

                tempMode_2 = reshape(shiftdim(tempCore,1),M(2),I(3)*I(1));
                tempResult = D{2}*tempMode_2;
                tempCore = shiftdim(reshape(tempResult,I(2),I(3),I(1)),2);

                Error = Y - tensor(tempCore);
                P = parafac_als(Error,1);
                D{3}(:,i) = P.U{3};                
            end
        end

        clear Xmode_3 tempXmode_3 tempMode_1 tempMode_2 Error tempResult tempCore;

        X=sptensor(X);
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
% æŠŠvä¸­çš„å�?nä¸ªæ•°å¾ªçŽ¯å�?å³ç§»åŠ¨1ä½ï¼Œå…¶ä»�?çš„N-nä¸ªæ•°æ�?°ä¿æŒä¸å˜ï¼�?
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