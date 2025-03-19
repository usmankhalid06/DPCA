function [srtd_Zt,srtd_Zs,ind_Zs]=sort_TSandSM_spatial(TC,SM,Zt,Zs,srcs) 
    ind_Zs = [0];
   for j=1:srcs
        [~, tmp_ind_Zs]  = sort(abs(corr(abs(SM(j,:)'),abs(Zs'))),'descend');
        if ismember(tmp_ind_Zs(1), ind_Zs)
            srtd_Zs(j,:) =  Zs(tmp_ind_Zs(2),:);
            srtd_Zt(:,j) = sign(corr(TC(:,j),Zt(:,tmp_ind_Zs(2))))*Zt(:,tmp_ind_Zs(2));  
            ind_Zs(j) = tmp_ind_Zs(2);
        else
            srtd_Zs(j,:) =  Zs(tmp_ind_Zs(1),:);
            srtd_Zt(:,j) = sign(corr(TC(:,j),Zt(:,tmp_ind_Zs(1))))*Zt(:,tmp_ind_Zs(1)); 
            ind_Zs(j) = tmp_ind_Zs(1);
        end
   end  

