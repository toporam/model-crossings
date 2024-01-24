function LLFfigs_PLUS(chooseDatab, chooseImCond, doFig5, doFig6, nFaces, BasisPath, Path2Save, barColor)

views = [0, 45, 90, 135, 180];

if strcmp(chooseImCond, 'pixel')
    for id = 1:nFaces
        for ori = 1:length(views)
            imageS = imread([BasisPath '/ImageDatasets/pixel/' chooseDatab 'Only' filesep chooseDatab num2str(views(ori)) '_' num2str(id) '.bmp']);
            imageLeft = imageS(:,1:round(end/2));
            mnImg(id,ori) = mean(imageS(:)); %#ok<*AGROW>
            LmnImg(id,ori) = mean(imageLeft(:));
            vrImg(id,ori) = var(double(imageS(:)));
            LvrImg(id,ori) = var(double(imageLeft(:)));
        end
    end
elseif strcmp(chooseImCond, 's1')
    load([BasisPath '/ImageDatasets/s1/s1Img_' chooseDatab 'Only.mat'], 'S1Img');
    for id = 1:nFaces
        for ori = 1:length(views)
            imageS = cell2mat(S1Img(id,ori));
            imageLeft = imageS(:,1:round(end/2));
            mnImg(id,ori) = mean(imageS(:));
            LmnImg(id,ori) = mean(imageLeft(:));
            vrImg(id,ori) = var(imageS(:));
            LvrImg(id,ori) = var(imageLeft(:));
        end
    end
end

for id = 1:nFaces
    vectMn = mnImg(id,:);
    betasMn(id,:) = glmQuadLin(vectMn);
    LvectMn = LmnImg(id,:);
    LbetasMn(id,:) = glmQuadLin(LvectMn);
    vectVr = vrImg(id,:);
    betasVr(id,:) = glmQuadLin(vectVr);
    LvectVr = LvrImg(id,:);
    LbetasVr(id,:) = glmQuadLin(LvectVr);
end

%% ADD BETAS & R2 HERE

% Get R2's from GLMx0x2 & GLMx0x4. The latter has R2s partitioned into
% a symm & and an antisymm componenet Nov 5, 2023. This should allow me to
% get the pending stats at the image level without the need to use the hack

for id = 1:nFaces
    vectMn_NEW = mnImg(id,:);
    [betasMn_NEW(id,:), R2Mn_NEW(id,:)] = glm_x0x2(vectMn_NEW);
    [betasMn_x0x4(id,:), R2Mn_x0x4(id,:)] = glm_x0x4(vectMn_NEW);

    LvectMn_NEW = LmnImg(id,:);
    [LbetasMn_NEW(id,:), LR2Mn_NEW(id,:)] = glm_x0x2(LvectMn_NEW);
    [LbetasMn_x0x4(id,:), LR2Mn_x0x4(id,:)] = glm_x0x4(LvectMn_NEW);

    vectVr_NEW = vrImg(id,:);
    [betasVr_NEW(id,:), R2Vr_NEW(id,:)] = glm_x0x2(vectVr_NEW);
    [betasVr_x0x4(id,:), R2Vr_x0x4(id,:)] = glm_x0x4(vectVr_NEW);

    LvectVr_NEW = LvrImg(id,:);
    [LbetasVr_NEW(id,:), LR2Vr_NEW(id,:)] = glm_x0x2(LvectVr_NEW);
    [LbetasVr_x0x4(id,:), LR2Vr_x0x4(id,:)] = glm_x0x4(LvectVr_NEW);
    
end

% check here 1) that betas are the same, that R^2 values make sense
% keyboard
%%

if doFig5
    mnImgMed = median(mnImg);
    upper= prctile(mnImg,75) - mnImgMed;
    lower= mnImgMed - prctile(mnImg,25);

    fig = figure;

    subplot(1,2,1), bar(mnImgMed, 'FaceColor',[barColor barColor barColor]), set(gca,'xticklabel',{'-90','-45','0','45','90'}), title([chooseDatab ' ' chooseImCond]), xlabel('Face orientation'), ylabel('Luminance (mean)');
    hold on
    er = errorbar((1:length(views)),mnImgMed,lower,upper);
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    axis square

    %------------------------------------------
    % fit & plot polynomial
    %------------------------------------------
    xX = 1:length(views); orderPoly=2;
    p = polyfit(xX,mnImgMed,orderPoly); % fit polynomial of second degree
    x2 = 1:.1:length(views);
    y2 = polyval(p,x2);
    plot(x2,y2,'LineWidth',2,'Color',[1, 0, 0, 0.5]);
    hold off
    %------------------------------------------

    vrImgMed = median(vrImg);
    upper= prctile(vrImg,75) - vrImgMed;
    lower= vrImgMed - prctile(vrImg,25);

    subplot(1,2,2), bar(vrImgMed, 'FaceColor',[barColor barColor barColor]), set(gca,'xticklabel',{'-90','-45','0','45','90'}), title([chooseDatab ' ' chooseImCond]), xlabel('Face orientation'), ylabel('Contrast (variance)');
    hold on
    er = errorbar((1:length(views)),vrImgMed,lower,upper);
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    axis square
    %------------------------------------------
    % fit & plot polynomial
    %------------------------------------------
    p = polyfit(xX,vrImgMed,orderPoly); % fit polynomial of second degree
    x2 = 1:.1:length(views);
    y2 = polyval(p,x2);
    plot(x2,y2,'LineWidth',2,'Color',[1, 0, 0, 0.5]);
    hold off
    %------------------------------------------

    print(fig, [Path2Save, 'Figure5'], '-depsc', '-r300')
    if ~exist([BasisPath '/DATA/Figure5']) %#ok<EXIST>
        mkdir([BasisPath '/DATA/Figure5'])
    end
    cd ([BasisPath '/DATA/Figure5'])
    save([chooseDatab '_' chooseImCond '_Mean.mat'], 'mnImg')
    save([chooseDatab '_' chooseImCond '_Variance.mat'], 'vrImg')
    save(['betas' chooseDatab '_' chooseImCond '_Mean.mat'], 'betasMn', ...
        'betasMn_NEW', 'R2Mn_NEW','betasMn_x0x4', 'R2Mn_x0x4') % Nov 3, 2023. For stats R^2 w/o hack (cf. mainComputeStats.m)
    save(['betas' chooseDatab '_' chooseImCond '_Variance.mat'], 'betasVr', ...
        'betasVr_NEW', 'R2Vr_NEW', 'betasVr_x0x4', 'R2Vr_x0x4') % Nov 3, 2023. For stats R^2 w/o hack (cf. mainComputeStats.m)

elseif doFig6
    LmnImgMed = median(LmnImg);
    upperL= prctile(LmnImg,75) - LmnImgMed;
    lowerL= LmnImgMed - prctile(LmnImg,25);

    fig = figure;

    subplot(1,2,1), bar(LmnImgMed, 'FaceColor',[barColor barColor barColor]), set(gca,'xticklabel',{'-90','-45','0','45','90'}), title(['Left half: ' chooseDatab ' ' chooseImCond]), xlabel('Face orientation'), ylabel('Luminance (mean)');
    hold on
    er = errorbar((1:length(views)),LmnImgMed,lowerL,upperL);
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    axis square
    %------------------------------------------
    % fit & plot polynomial
    %------------------------------------------
    xX = 1:length(views); orderPoly=2;
    p = polyfit(xX,LmnImgMed,orderPoly); % fit polynomial of second degree
    x2 = 1:.1:length(views);
    y2 = polyval(p,x2);
    plot(x2,y2,'LineWidth',2,'Color',[1, 0, 0, 0.5]);
    hold off
    %------------------------------------------

    LvrImgMed = median(LvrImg);
    upperL= prctile(LvrImg,75) - LvrImgMed;
    lowerL= LvrImgMed - prctile(LvrImg,25);

    subplot(1,2,2), bar(LvrImgMed, 'FaceColor',[barColor barColor barColor]), set(gca,'xticklabel',{'-90','-45','0','45','90'}), title(['Left half: ' chooseDatab ' ' chooseImCond]), xlabel('Face orientation'), ylabel('Contrast (variance)');
    hold on
    er = errorbar((1:length(views)),LvrImgMed,lowerL,upperL);
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    axis square
    %------------------------------------------
    % fit & plot polynomial
    %------------------------------------------
    p = polyfit(xX,LvrImgMed,orderPoly); % fit polynomial of second degree
    x2 = 1:.1:length(views);
    y2 = polyval(p,x2);
    plot(x2,y2,'LineWidth',2,'Color',[1, 0, 0, 0.5]);
    hold off
    %------------------------------------------

    hold off

    print(fig, [Path2Save, 'Figure6'], '-depsc', '-r300')
    if ~exist([BasisPath '/DATA/Figure6']) %#ok<EXIST>
        mkdir([BasisPath '/DATA/Figure6'])
    end
    cd ([BasisPath '/DATA/Figure6'])
    save([chooseDatab '_' chooseImCond '_Mean.mat'], 'LmnImg')
    save([chooseDatab '_' chooseImCond '_Variance.mat'], 'LvrImg')
    save(['betas' chooseDatab '_' chooseImCond '_Mean.mat'], 'LbetasMn', ...
        'LbetasMn_NEW', 'LR2Mn_NEW', 'LbetasMn_x0x4', 'LR2Mn_x0x4') % Nov 3, 2023. For stats R^2 w/o hack (cf. mainComputeStats.m))
    save(['betas' chooseDatab '_' chooseImCond '_Variance.mat'], 'LbetasVr', ...
        'LbetasVr_NEW','LR2Vr_NEW', 'LbetasVr_x0x4','LR2Vr_x0x4') % Nov 3, 2023. For stats R^2 w/o hack (cf. mainComputeStats.m))
end

%% add required subfunction
    function [betasAll, R2All] = glm_x0x2(vects)

        x = [-2 -1 0 1 2];
        y = x-mean(x);
        r1 = y/std(y);

        y = x.^2;
        y = y - mean(y);
        r2 = y/std(y);

        X = [r1', r2'];

        betasAll = glmfit(X,vects');

        symmComp = X*([0, betasAll(3)])';
        varSymmComp = var(symmComp);

        antisymmComp = X*([betasAll(2), 0])';
        varAntisymmComp = var(antisymmComp);

        totalVar = var(vects);

        varSigRecon = var(X*([betasAll(2), betasAll(3)])');

        checkSumVarRecon = varSigRecon-(varAntisymmComp+varSymmComp); %#ok<*NASGU>

        Res = vects-(X*([betasAll(2), betasAll(3)])')';
        varRes = var(Res);

        R2All = [varAntisymmComp/totalVar, varSymmComp/totalVar];

        checkSumR2 = sum(R2All)+(varRes/totalVar);

    end

    function [betasAll, R2All] = glm_x0x4(vects) %#ok<*DEFNU>
        % GLM with antisymmetric (linear + cubic) and symmetric (quadratic +
        % quartic) regressors

        x = [-2 -1 0 1 2];
        y = x-mean(x);
        r1 = y/std(y);

        y = x.^2;
        y = y - mean(y);
        r2 = y/std(y);

        y = x.^3;
        y = y - mean(y);
        r3 = y/std(y);

        y = x.^4;
        y = y - mean(y);
        r4 = y/std(y);

        X = [r1', r2', r3', r4'];

        betasAll = glmfit(X,vects');

        symmComp = X*([0, betasAll(3), 0, betasAll(5)])';
        varSymmComp = var(symmComp);

        antisymmComp = X*([betasAll(2), 0, betasAll(4), 0])';
        varAntisymmComp = var(antisymmComp);

        totalVar = var(vects);

        R2All = [varAntisymmComp/totalVar, varSymmComp/totalVar];
        checkSumR2 = sum(R2All);

    end

end

