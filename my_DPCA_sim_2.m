function [D,X,Err,CC,A,B]= my_DPCA_sim_2(Y,Dp,Xp,K,spa1,nIter,TC,SM)
    K1 = size(Dp,2);
    K2 = size(Xp,1);
    A = eye(K1,K);
    B = zeros(K,K2);
    D = Dp*A;
    X = B*Xp;
    sizeD = numel(D);   
    D_tol = 1e-9;

%     param.lambda=spa1;
%     param.numThreads=4; % number of threads
%     param.iter=nIter;     
    %     fprintf('Iteration:     ');
    for iter=1:nIter
        %         fprintf('\b\b\b\b\b%5i',iter);
        Do = D;

%         xkk=full(mexLasso(Y,D,param));
%         B = xkk*pinv(Xp);
%         X = B*Xp;
       
            for i =1:K
                xk = D(:,i)'*Y;
                thr1 = spa1./abs(xk); 
                xkk = sign(xk).*max(0, bsxfun(@minus,abs(xk),thr1/2));
                B(i,:) = xkk*pinv(Xp);
                X(i,:) = B(i,:)*Xp;
            end     
%         
        F2 = X*X'; E2 = Y*X';
        iiter = 0;
        Dpp = D; 
        while (iiter < nIter)
            iiter = iiter + 1;
            for j = 1: K
                if(F2(j,j) ~= 0)
                    tmp3 = 1.0/F2(j,j) * (E2(:,j) - D*F2(:, j)) + D(:,j);
                    A(:,j)= pinv(Dp)*tmp3; 
%                     A(:,j)  = TLS(Dp, tmp3); 
                    A(:,j) = A(:,j)./(max( norm(Dp*A(:,j),2),1));
                    D(:,j) = Dp*A(:,j);
                end
            end
            %% check stop condition
            if (norm(D - Dpp, 'fro')/sizeD < D_tol)
                break;
            end
            Dpp = D;
        end
        Err(iter)= sqrt(trace((D-Do)'*(D-Do)))/sqrt(trace(Do'*Do));
    
    
    
        [~,~,ind]=sort_TSandSM_spatial(TC,SM,D,X,K);
        for ii =1:K
            TCcorr(ii) =abs(corr(TC(:,ii),D(:,ind(ii))));
            SMcorr(ii) =abs(corr(SM(ii,:)',X(ind(ii),:)'));
        end
        cTC = sum(TCcorr');
        cSM = sum(SMcorr');
        CC(iter) =cTC+cSM;
    end

%     if ref ==1
%         spa =spa/2;
%         X= sign(X).*max(0, bsxfun(@minus,abs(X),spa/2));
%         D = Y*pinv(X);
%         D = D*diag(1./sqrt(sum(D.*D))); 
%     end

end
