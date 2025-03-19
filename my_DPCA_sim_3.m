function [D,X,Err,A,B]= my_dPCA_sim_3(Y,Dp,Xp,K,spa1,spa2,spa3,nIter)
    K1 = size(Dp,2);
    K2 = size(Xp,1);
    A = eye(K1,K);
    B = zeros(K,K2);
    D = Dp*A;
    X = B*Xp;
    fprintf('Iteration:     ');
    for iter=1:nIter
        fprintf('\b\b\b\b\b%5i',iter);
        Do = D;
        for j =1:K
            X(j,:) = 0; A(:,j) = 0; B(j,:) = 0;
            E = Y-D*X;
            
            xk = D(:,j)'*E;
            thr1 = spa1./abs(xk); 
            xkk = sign(xk).*max(0, bsxfun(@minus,abs(xk),thr1/2));
%             B(j,:) = TLS(Xp', xkk');
            B(j,:) = xkk*pinv(Xp);
            X(j,:) = B(j,:)*Xp;

%             X(j,:) = sign(X(j,:)).*max(0, bsxfun(@minus,abs(X(j,:)),thr1/2));
            X(j,:) = firm_thresholding_nonadaptive(X(j,:), spa2/2, spa3/2);  
%             B(j,:) = X(j,:)*pinv(Xp);
%             X(j,:) = B(j,:)*Xp;            
            
            rInd = find(X(j,:));
            if (length(rInd)>1)
                tmp3 = E(:,rInd)*X(j,rInd)'; 
                A(:,j)= pinv(Dp)*tmp3;    
%                 A(:,j)  = TLS(Dp, tmp3); 
                A(:,j) = A(:,j)./norm(Dp*A(:,j));
                D(:,j) = Dp*A(:,j);
            end          
        end     
        Err(iter) = sqrt(trace((D-Do)'*(D-Do)))/sqrt(trace(Do'*Do));

    end
end


