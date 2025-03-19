clear;
close all; 
clc;



%% parameters
sp = simtb_create_sP('exp_params_aod');
N = sp.nT; %number of time samples 
nV = sp.nV; %sqrt of number of voxels
nSRCS = sp.nC; %number of sources
K =8 ; %dimensionality reduction sources
nIter = 30; %algorithm iterations
tstd  = sqrt(0.9); %0.9
sstd  =  sqrt(0.005); %0.01
Dp = dctbases(N,N); %dct basis dictionary
nf = sqrt(70*15.49); 

% for j=1:1
j = 1;    

    sp.SM_spread =8+0.05*randn(N,nV*nV); %0.0001
    SM = simtb_makeSM(sp,1);   % Create spatial maps
    TC = zscore(simtb_makeTC(sp,1));  % Create TCs 

    
    TC = zscore(Dp(:,3:8:60));
    rng('default'); 
    rng('shuffle') % random number generator
    Y= (TC+tstd*randn(N,nSRCS))*(SM+sstd*randn(nSRCS,nV^2));
    Y= Y-repmat(mean(Y),size(Y,1),1);    

    %% sICA
    tStart=tic;
    [G,~,~] = svds(Y,K);
    Ss = G'*Y;
    [SSs,A,~] = fastica(Ss, 'numOfIC', K-0,'approach','symm', 'g', 'tanh','verbose', 'off');
    X{1} = SSs;
    D{1} = Y*SSs'; 
    D{1} = D{1}*diag(1./sqrt(sum(D{1}.*D{1})));
    tEnd(j,1) = toc(tStart);  

    %% Karim SPCA
    tStart=tic;
    [~,~,~,Ss] = Sparse_PCA(Y, 0.005, K, nIter, nIter, 1e-2, 1e-2);
    [SSs,A,~] = fastica(Ss', 'numOfIC', K-0,'approach','symm', 'g', 'tanh','verbose', 'off');
    X{2} = SSs;
    D{2} = Y*SSs'; 
    D{2} = D{2}*diag(1./sqrt(sum(D{2}.*D{2})));
    tEnd(j,2) = toc(tStart);      

    %% SPCA
    tStart=tic;
    for i =1:K;  lambda(i)=sum((SM(i,:)>=1e-1));  end
    [tmpX,~] = spca(Y,[],K,Inf,-lambda); % 2000
    [SSs,A,~] = fastica(tmpX', 'numOfIC', K-0,'approach','symm', 'g', 'tanh','verbose', 'off');
    X{3} = SSs;
    D{3} = Y*SSs';
    D{3} = D{3}*diag(1./sqrt(sum(D{3}.*D{3})));
    tEnd(j,3) = toc(tStart);
    
    %% SPC
    tStart=tic;
    lambda = 0.43; %0.43*70=31.5
    [Wx,~,di] = pmd_rankK(Y',Y',K,lambda);
    X{4} = Wx';
    D{4} = Y*X{4}';
    D{4} = D{4}*diag(1./sqrt(sum(D{4}.*D{4})));
    tEnd(j,4) = toc(tStart);

    %% GPower
    tStart=tic;
    lambda =0.3;
    gamma=lambda*ones(1,K);
    Wx=GPower(Y,gamma,K,'l1',0);
    X{5} = Wx';
    D{5} = Y*X{5}';
    D{5} = D{5}*diag(1./sqrt(sum(D{5}.*D{5})));
    tEnd(j,5) = toc(tStart);

    %% DPCA2
    tStart=tic;
    spa1 = (22/nf)*nf; %18 1
    [Dq, V, Xq]= svds(Y,K);    
    [D{6},X{6},Err,CC,A2,B2]= my_DPCA_sim_2(Y,Dq,Xq',K,spa1,nIter,TC,SM);
    tEnd(j,6) = toc(tStart);
%     figure; plot(Err)

    %% DPCA1
    tStart=tic;
    spa1 = (22/nf)*nf; 
    spa2 = (4/nf)*nf;  %6
    spa3 = (8/nf)*nf; %16
    [Dq, V, Xq]= svds(Y,K);    
    [D{7},X{7},Err,CC,A1,B1]= my_DPCA_sim_1(Y,Dq,Xq',K,spa1,spa2,spa3,nIter,TC,SM);
    tEnd(j,7) = toc(tStart);
%     figure; plot(Err)   

    %% DPCA3
    tStart=tic;
    spa1 =(22/nf)*nf; %40 (for no V included in V*Xq')  15 (otherwise)
    spa2 = (4/nf)*nf; %120 (for no V included in V*Xq')  60 (otherwise)
    spa3 = (8/nf)*nf;
    [Dq, V, Xq]= svds(Y,K); 
    [D{8},X{8},Err,A3,B3]= my_DPCA_sim_3(Y,Dq,Xq',K,spa1,spa2,spa3,nIter);
    tEnd(j,8) = toc(tStart);
%     figure; plot(Err)
   
    %%
    sD{1} = TC; 
    sX{1} = SM;
    nA = 9;
    j=1;
    for jj =1:nA-1
        [sD{jj+1},sX{jj+1},ind{jj}]=sort_TSandSM_temporal(TC,D{jj},X{jj});
%         [sD{jj+1},sX{jj+1},ind{jj}]=sort_TSandSM_spatial(TC,SM,D{jj},X{jj},nSRCS);
        for ii =1:nSRCS
            TCcorr(jj+1,ii,j) =abs(corr(TC(:,ii),D{jj}(:,ind{jj}(ii))));
            SMcorr(jj+1,ii,j) =abs(corr(SM(ii,:)',X{jj}(ind{jj}(ii),:)'));
        end
    end


% end


f = figure; f.Position = [10 50 650 950]; nA = 8;my_subplots_horz2(nA+1,K,nV,nV,TCcorr,SMcorr,sD,sX);
exportgraphics(gcf,'fig1b.png','Resolution',300)

