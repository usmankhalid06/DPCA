function my_subplots_horz2(nA,S,v,w,TCcorr,SMcorr,rTC,rSM)

    axis off
    set(gca,'Units','normalized','Position',[0 0 1 1]);
    str_title_v={'(A)','(B)','(C)','(D)','(E)','(F)','(G)','(H)','(I)'};
    str_title_h={'lv_1','lv_2','lv_3','lv_4','lv_5','lv_6','lv_7','lv_8'};


    a_vec = 0.93:-0.105:0.05;
    b_vec = 0.92:-0.115:0.01;
    text(0.01,a_vec(1), '(A)','Color','b','FontSize',10)
    text(b_vec(1), 0.98, str_title_h{1},'Color','b','FontSize',10)
    for i =2:9
        text(0.01,a_vec(i)-0.055, [num2str(round(mean(SMcorr(i,:)),2),'%0.2f') ],'Color','m','FontSize',10)
        text(0.01,a_vec(i), str_title_v{i},'Color','b','FontSize',10)
        if i<=8
         text(b_vec(i), 0.98, str_title_h{i},'Color','b','FontSize',10)
        end
    end

    centerX = v/2; % X-coordinate of the circle center (adjust as needed)
    centerY = v/2; % Y-coordinate of the circle center (adjust as needed)
    radius = v/2;
    [x, y] = meshgrid(1:w, 1:v);
    distanceFromCenter = sqrt((x - centerX).^2 + (y - centerY).^2);
    circleMask = distanceFromCenter <= radius;
       
    ivs = 1.34;  %initial_horizontal_shift
    ihs = 1.3;  %initial_vertical_shift
    shz = S-1.00; %subplot_horizontal_size
    svs = 0.1;  %subplot_vertical_size (more the better)
    vs  = 0.65;  %vertical shift of subplots
    nC  = S+0.75; %No. of columns
    shifter = 0.16; %vertical
    for i =1:nA 
        for j=1:S

            if i==1
            %%
                zscore_rxSM = abs(zscore(rSM{i}(j,:)));
                activationImage = flipdim(reshape(zscore_rxSM,v,w),1);
                maskedActivationImage = activationImage .* circleMask;
                maskedActivationImage(~circleMask) = 0;                

                hax=axes();
                set(gca, 'Color', 'k'); hold on;
                imagesc(maskedActivationImage);  
                cmap = jet(256);
                cmap(1,:) = 0;
                colormap(cmap);
                rectangle('Position', [centerX - radius, centerY - radius, 2 * radius, 2 * radius], ...
                'Curvature', [1, 1], 'EdgeColor', 'white', 'LineWidth', 1);

                newPos=[(1-1/nC)-(1/nC)*(fix((j-1)/1)+ihs-1),  vs*(mod(j-1,1)+ivs+0.00),     1/shz,   svs];
                set(gca,'outer',newPos), 
                set(gca,'XTickLabel','')
                set(gca,'YTickLabel','')
%                 
            else

            %%
                zscore_rxSM = abs(zscore(rSM{i}(j,:)));
                activationImage = flipdim(reshape(zscore_rxSM,v,w),1);
%                 activationImage(activationImage<=0.001) = 100;
                maskedActivationImage = activationImage .* circleMask;
                maskedActivationImage(~circleMask) = 0;

                hax=axes();
                set(gca, 'Color', 'k'); hold on;
                imagesc(maskedActivationImage);  
                cmap = jet(256);
                cmap(1,:) = 0;
                colormap(cmap);
                rectangle('Position', [centerX - radius, centerY - radius, 2 * radius, 2 * radius], ...
                'Curvature', [1, 1], 'EdgeColor', 'white', 'LineWidth', 1);

                newPos=[(1-1/nC)-(1/nC)*(fix((j-1)/1)+ihs-1),  vs*(mod(j-1,1)+ivs-shifter*(i-1)),     1/shz,   svs];
                set(gca,'outer',newPos), 
                set(gca,'XTickLabel','')
                set(gca,'YTickLabel','')
                xlabel([num2str(round(SMcorr(i,j),2))],'color','k')
                xh = get(gca,'xlabel'); % handle to the label object
                p = get(xh,'position'); % get the current position property
                p(2) = 5+ p(2);       % double the distance,
                set(xh,'position',p)   % set the new position    

            end
        end
    end
