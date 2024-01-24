function callCrossingsModel(MODEL, BasisPath)

rngSeed = rng('default'); %#ok<NASGU>
scrsz = get(groot, 'Screensize');

%-------------------------------------------------------------------------
% Define grids for generation of probability distributions
%-------------------------------------------------------------------------
start = -6; step = 12/440; stop = 6; % parameters for mesh grid
vMeshGrid = [start, step, stop];
%-------------------------------------------------------------------------
% Prepare auxiliary vectors and variables used in script
%-------------------------------------------------------------------------
vDatab = {'radboud','kdef'};
vImCond = {'pixel','s1'};
vCortMag = {'CM','noCM'};

chooseDatab = vDatab{MODEL.indDatab};
chooseImCond = vImCond{MODEL.indImCond};
chooseCortMag = vCortMag{MODEL.indCortMag};
nFaces = size(dir([BasisPath '/ImageDatasets/pixel/' chooseDatab filesep chooseDatab '0*.bmp']),1); %number of IDs in the database

if MODEL.applyGF == 0
    chooseGF = 'noGF';
elseif MODEL.applyGF == 1
    chooseGF = 'GF';
end

%-------------------------------------------------------------------------
% Define path2save
%-------------------------------------------------------------------------
JobName = [chooseDatab '_' chooseImCond '_' chooseCortMag '_' chooseGF];

Path2Save = [BasisPath, filesep, 'Results', filesep, JobName, filesep];

if ~exist(Path2Save)
    mkdir(Path2Save)
end

%% FIGURE 4A

[Z4R, Z4Rvect, Z4L, Z4Lvect] = mkDists(MODEL.k, vMeshGrid, MODEL.doFig4, MODEL.indCortMag, BasisPath,scrsz);

if MODEL.doFig4
    return
end

%% Prepare images

if MODEL.indImCond == 2 % if S1 model variant
    if ~exist([BasisPath '/ImageDatasets/s1'])
        cd([BasisPath '/ImageDatasets/'])
        mkdir(chooseImCond)
    end

    if ~exist([BasisPath '/ImageDatasets/s1/s1Img_' chooseDatab '.mat'])
        imgAsIs = 0;
        S1Img = getS1(chooseDatab, nFaces, BasisPath, MODEL.seeMadeImg, imgAsIs);
        cd ([BasisPath '/ImageDatasets/s1/'])
        save(['s1Img_' chooseDatab '.mat'],'S1Img')
    end

    if ~exist([BasisPath '/ImageDatasets/s1/s1Img_' chooseDatab 'Only.mat'])
        imgAsIs = 1;
        S1Img = getS1(chooseDatab, nFaces, BasisPath, MODEL.seeMadeImg, imgAsIs);
        cd ([BasisPath '/ImageDatasets/s1/'])
        save(['s1Img_' chooseDatab 'Only.mat'],'S1Img')
    end
end

%% FIGURE 5 AND 6

if MODEL.doFig5 || MODEL.doFig6
    % LLFfigs(chooseDatab, chooseImCond, MODEL.doFig4, MODEL.doFig5, nFaces, BasisPath, Path2Save, MODEL.barColor);
    LLFfigs_PLUS(chooseDatab, chooseImCond, MODEL.doFig5, MODEL.doFig6, nFaces, BasisPath, Path2Save, MODEL.barColor);    
    return
end

%% Specify network connections, or load pre-specified network

if MODEL.mkConnect
    indxL1 = mkL1(Z4Rvect, Z4Lvect, MODEL.nUnitsL1, MODEL.density); % layer 1
    indxLRest = mkLRest(MODEL.nUnitsL1, MODEL.nUnitsRest, MODEL.density, MODEL.nLayers, MODEL.biProb); % layers 2-8
    cd([BasisPath '/Results'])
    save(['netSettings_' chooseCortMag '.mat'],'indxL1','indxLRest')
else
    cd([BasisPath '/Results'])
    load(['netSettings_' chooseCortMag '.mat'],'indxL1','indxLRest')
end

%% Generate gain field

if ~exist([BasisPath '/Results/gainfield.mat']) %#ok<*EXIST>
    GFL1 = rand(MODEL.nUnitsL1*2,1); 
    allGF{1} = repmat(GFL1,1,5,length(MODEL.density));

    GFLrest = rand(MODEL.nUnitsRest,MODEL.nLayers-1);
    for i = 1:MODEL.nLayers-1
        allGF{i+1} = repmat(GFLrest(:,i),1,5,length(MODEL.density)); 
    end
    cd([BasisPath '/Results'])
    save('gainfield.mat','allGF')
else
    cd([BasisPath '/Results'])
    load('gainfield.mat', 'allGF')
end

%% Get activation patterns from pre-specified network

if ~exist([Path2Save '/ActivPattern.mat'])

    for id = 1:nFaces %go through all IDs
        [imageV] = prepImgEvent(chooseDatab, chooseImCond, id, MODEL.seeImg, BasisPath);

        [L1Actis, RestActisC] = defineActiv(imageV, indxL1, indxLRest, MODEL.nLayers, MODEL.density);

        allLayers{1} = L1Actis;
        allLayers{2} = RestActisC{1};
        allLayers{3} = RestActisC{2};
        allLayers{4} = RestActisC{3};
        allLayers{5} = RestActisC{4};
        allLayers{6} = RestActisC{5};
        allLayers{7} = RestActisC{6};
        allLayers{8} = RestActisC{7};

        % apply GF
        if MODEL.applyGF == 1
            allLayers = applyGF(allLayers,allGF);
        end

        PatternSaved{id,:} = allLayers; %save patterns for all images
        disp(['ID ' num2str(id) ' done'])

    end

    cd (Path2Save)
    save('ActivPattern.mat', 'PatternSaved')
else
    cd (Path2Save)
    load('ActivPattern.mat', 'PatternSaved')
end

%% FIGURE 7

if MODEL.doFig7

    if ~exist([BasisPath '/DATA/Figure7/' chooseCortMag filesep chooseGF filesep chooseDatab '_' chooseImCond filesep 'betasR2.mat'])

        for lay = 1:MODEL.nLayers

            for id = 1:nFaces

                for dens = 1:length(MODEL.density)
                    aux = PatternSaved{id}{lay};
                    aux = aux(:,:,dens);
                    halfUnits = length(aux)/2;

                    for ori = 1:5
                        vectMn(ori) = mean(aux(:,ori)); %#ok<*AGROW> %full layer
                        vectVr(ori) = var(aux(:,ori));

                        vectMnLeft(ori) = mean(aux(1:halfUnits,ori)); %left half of layer
                        vectVrLeft(ori) = var(aux(1:halfUnits,ori));

                        vectMnRight(ori) = mean(aux(halfUnits+1:end,ori)); %right half of layer
                        vectVrRight(ori) = var(aux(halfUnits+1:end,ori));
                    end

                    meanLay(id,:,lay,dens)= vectMn;
                    varLay(id,:,lay,dens) = vectVr;

                    [betasMnLocal, r2MnLocal] = glm_x0x4(vectMn);
                    [betasVrLocal, r2VrLocal] = glm_x0x4(vectVr);

                    betasMn(:,lay,dens,id) = betasMnLocal;
                    r2Mn(:,lay,dens,id) = r2MnLocal;
                    betasVr(:,lay,dens,id) = betasVrLocal;
                    r2Vr(:,lay,dens,id) = r2VrLocal;

                    meanLeftLay(id,:,lay,dens)= vectMnLeft; %#ok<NASGU>
                    varLeftLay(id,:,lay,dens) = vectVrLeft; %#ok<NASGU>

                    [betasMnLeftLocal, r2MnLeftLocal] = glm_x0x4(vectMnLeft);
                    [betasVrLeftLocal, r2VrLeftLocal] = glm_x0x4(vectVrLeft);

                    betasMnLeft(:,lay,dens,id) = betasMnLeftLocal;
                    r2MnLeft(:,lay,dens,id) = r2MnLeftLocal;
                    betasVrLeft(:,lay,dens,id) = betasVrLeftLocal;
                    r2VrLeft(:,lay,dens,id) = r2VrLeftLocal;

                    meanRightLay(id,:,lay,dens)= vectMnRight;
                    varRightLay(id,:,lay,dens) = vectVrRight;

                    [betasMnRightLocal, r2MnRightLocal] = glm_x0x4(vectMnRight);
                    [betasVrRightLocal, r2VrRightLocal] = glm_x0x4(vectVrRight);

                    betasMnRight(:,lay,dens,id) = betasMnRightLocal;
                    r2MnRight(:,lay,dens,id) = r2MnRightLocal;
                    betasVrRight(:,lay,dens,id) = betasVrRightLocal;
                    r2VrRight(:,lay,dens,id) = r2VrRightLocal;

                end
            end
        end

        %----------------------------------------------------------------------
        if ~exist([BasisPath '/DATA/Figure7/' chooseCortMag filesep chooseGF filesep chooseDatab '_' chooseImCond])
            mkdir([BasisPath '/DATA/Figure7/' chooseCortMag filesep chooseGF filesep chooseDatab '_' chooseImCond])
        end

        cd([BasisPath '/DATA/Figure7/' chooseCortMag filesep chooseGF filesep chooseDatab '_' chooseImCond])

        if MODEL.regressFit == 1 % difference of betas
            for dens = 1:length(MODEL.density)
                for id = 1:nFaces
                    diffMnFull(id,dens) = abs(betasMn(3,MODEL.pickLay,dens,id)) - abs(betasMn(2,MODEL.pickLay,dens,id));
                    diffMnLeft(id,dens) = abs(betasMnLeft(3,MODEL.pickLay,dens,id)) - abs(betasMnLeft(2,MODEL.pickLay,dens,id));
                    diffMnRight(id,dens) = abs(betasMnRight(3,MODEL.pickLay,dens,id)) - abs(betasMnRight(2,MODEL.pickLay,dens,id));
                    diffVrFull(id,dens) = abs(betasVr(3,MODEL.pickLay,dens,id)) - abs(betasVr(2,MODEL.pickLay,dens,id));
                    diffVrLeft(id,dens) = abs(betasVrLeft(3,MODEL.pickLay,dens,id)) - abs(betasVrLeft(2,MODEL.pickLay,dens,id));
                    diffVrRight(id,dens) = abs(betasVrRight(3,MODEL.pickLay,dens,id)) - abs(betasVrRight(2,MODEL.pickLay,dens,id));
                end
            end
            save('betadiffMean.mat', 'diffMnFull', 'diffMnLeft', 'diffMnRight')
            save('betadiffVariance.mat', 'diffVrFull', 'diffVrLeft', 'diffVrRight')

        elseif MODEL.regressFit == 2 % difference of R^2's

            for dens = 1:length(MODEL.density)
                for id = 1:nFaces

                    diffMnFull(id,dens) = r2Mn(2,MODEL.pickLay,dens,id) - r2Mn(1,MODEL.pickLay,dens,id);
                    diffMnLeft(id,dens) = r2MnLeft(2,MODEL.pickLay,dens,id) - r2MnLeft(1,MODEL.pickLay,dens,id);
                    diffMnRight(id,dens) = r2MnRight(2,MODEL.pickLay,dens,id) - r2MnRight(1,MODEL.pickLay,dens,id);
                    diffVrFull(id,dens) = r2Vr(2,MODEL.pickLay,dens,id) - r2Vr(1,MODEL.pickLay,dens,id);
                    diffVrLeft(id,dens) = r2VrLeft(2,MODEL.pickLay,dens,id) - r2VrLeft(1,MODEL.pickLay,dens,id);
                    diffVrRight(id,dens) = r2VrRight(2,MODEL.pickLay,dens,id) - r2VrRight(1,MODEL.pickLay,dens,id);
                end
            end

            save('r2diffMean.mat', 'diffMnFull', 'diffMnLeft', 'diffMnRight')
            save('r2diffVariance.mat', 'diffVrFull', 'diffVrLeft', 'diffVrRight')
        end

        save('betasR2.mat', 'meanLay', 'varLay', 'meanRightLay', 'varRightLay')

    else
        cd ([BasisPath '/DATA/Figure7/' chooseCortMag filesep chooseGF filesep chooseDatab '_' chooseImCond])
        if MODEL.regressFit == 1
            load('betadiffMean.mat') %#ok<*LOAD>
            load('betadiffVariance.mat')
        elseif MODEL.regressFit == 2
            load('r2diffMean.mat')
            load('r2diffVariance.mat')
        end
        load('betasR2.mat')
    end

    % panel a/b
    fig = figure;

    meanLayAux = meanLay(:,:, MODEL.pickLay,MODEL.pickDens);
    meanLayMed = median(meanLayAux);
    upper= prctile(meanLayAux,75) - meanLayMed;
    lower= meanLayMed - prctile(meanLayAux,25);

    subplot(2,2,1), bar(meanLayMed, 'FaceColor',[.9 .9 .9]), set(gca,'xticklabel',{'-90','-45','0','45','90'}), title(['Full layer: ' chooseDatab ' ' chooseImCond ' ' chooseCortMag]), xlabel('Face orientation'), ylabel('Pattern activation (mean n.u.)')
    hold on
    er = errorbar((1:5),meanLayMed,lower,upper);
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    axis square

    %------------------------------------------
    % fit & plot polynomial
    %------------------------------------------
    x = 1:5; orderPoly=2;
    p = polyfit(x,meanLayMed,orderPoly); % fit polynomial of second degree
    x2 = 1:.1:5;
    y2 = polyval(p,x2);
    plot(x2,y2,'LineWidth',2,'Color',[1, 0, 0, 0.5]);
    hold off
    %------------------------------------------

    varLayAux = varLay(:,:, MODEL.pickLay,MODEL.pickDens);
    varLayMed = median(varLayAux);
    upper= prctile(varLayAux,75) - varLayMed;
    lower= varLayMed - prctile(varLayAux,25);

%     subplot(2,2,2), bar(varLayMed, 'FaceColor',[.1 .1 .1]),set(gca,'xticklabel',{'-90','-45','0','45','90'}), title(['Full layer: ' chooseDatab ' ' chooseImCond ' ' chooseCortMag]), xlabel('Face orientation'), ylabel(["Pattern contrast (variance of n.u.'s)")
    subplot(2,2,2), bar(varLayMed, 'FaceColor',[.1 .1 .1]),set(gca,'xticklabel',{'-90','-45','0','45','90'}), title(['Full layer: ' chooseDatab ' ' chooseImCond ' ' chooseCortMag]), xlabel('Face orientation'), ylabel({'Pattern contrast (variance of n.u.''s)'})
    hold on
    er = errorbar((1:5),varLayMed,lower,upper);
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    axis square

    %------------------------------------------
    % fit & plot polynomial
    %------------------------------------------
    x = 1:5; orderPoly=2;
    p = polyfit(x,varLayMed,orderPoly); % fit polynomial of second degree
    x2 = 1:.1:5;
    y2 = polyval(p,x2);
    plot(x2,y2,'LineWidth',2,'Color',[1, 0, 0, 0.5]);
    hold off
    %------------------------------------------

    meanRightLayAux = meanRightLay(:,:, MODEL.pickLay,MODEL.pickDens);
    meanRightLayMed = median(meanRightLayAux);
    upper= prctile(meanRightLayAux,75) - meanRightLayMed;
    lower= meanRightLayMed - prctile(meanRightLayAux,25);

    subplot(2,2,3), bar(meanRightLayMed, 'FaceColor',[.9 .9 .9]), set(gca,'xticklabel',{'-90','-45','0','45','90'}), title(['Right half layer: ' chooseDatab ' ' chooseImCond ' ' chooseCortMag]), xlabel('Face orientation'), ylabel('Pattern activation (mean n.u.)')
    hold on
    er = errorbar((1:5),meanRightLayMed,lower,upper);
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    axis square

    %------------------------------------------
    % fit & plot polynomial
    %------------------------------------------
    x = 1:5; orderPoly=2;
    p = polyfit(x,meanRightLayMed,orderPoly); % fit polynomial of second degree
    x2 = 1:.1:5;
    y2 = polyval(p,x2);
    plot(x2,y2,'LineWidth',2,'Color',[1, 0, 0, 0.5]);
    hold off
    %------------------------------------------

    varRightLayAux = varRightLay(:,:, MODEL.pickLay,MODEL.pickDens);
    varRightLayMed = median(varRightLayAux);
    upper= prctile(varRightLayAux,75) - varRightLayMed;
    lower= varRightLayMed - prctile(varRightLayAux,25);

%     subplot(2,2,4), bar(varRightLayMed, 'FaceColor',[.1 .1 .1]), set(gca,'xticklabel',{'-90','-45','0','45','90'}), title(['Right half layer: ' chooseDatab ' ' chooseImCond ' ' chooseCortMag]), xlabel('Face orientation'), ylabel("Pattern contrast (variance of n.u.'s)")
    subplot(2,2,4), bar(varRightLayMed, 'FaceColor',[.1 .1 .1]), set(gca,'xticklabel',{'-90','-45','0','45','90'}), title(['Right half layer: ' chooseDatab ' ' chooseImCond ' ' chooseCortMag]), xlabel('Face orientation'), ylabel({'Pattern contrast (variance of n.u.''s)'})
    hold on
    er = errorbar((1:5),varRightLayMed,lower,upper);
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    axis square

    %------------------------------------------
    % fit & plot polynomial
    %------------------------------------------
    x = 1:5; orderPoly=2;
    p = polyfit(x,varRightLayMed,orderPoly); % fit polynomial of second degree
    x2 = 1:.1:5;
    y2 = polyval(p,x2);
    plot(x2,y2,'LineWidth',2,'Color',[1, 0, 0, 0.5]);
    hold off
    %------------------------------------------

    print(fig, [Path2Save, 'Figure7ab'], '-depsc', '-r300')

    if MODEL.indCortMag == 1
        return
    end

    %panel c
    if MODEL.regressFit == 1
        regressNm = 'beta';
    elseif MODEL.regressFit == 2
        regressNm = 'R^2';
    end

    fig = figure;

    Xaxis = 1:length(MODEL.density);
    s = [Xaxis, fliplr(Xaxis)];

    cd([BasisPath '/DATA/Figure7/CM' filesep chooseGF filesep chooseDatab '_' chooseImCond filesep])
    if MODEL.regressFit == 1
        load('betadiffMean.mat')
        load('betadiffVariance.mat')
    elseif MODEL.regressFit == 2
        if exist('r2diffMean.mat') && exist('r2diffVariance.mat')
            load('r2diffMean.mat')
            load('r2diffVariance.mat')
        else
            error('NON EXISTING MAT-FILE(S) r2diffVariance.mat OR r2diffMean.mat not found')
        end
    end

    ci = bootci(2000,@mean,diffMnFull);
    lower = ci(1,:); upper = ci(2,:);
    btw = [lower, fliplr(upper)];
    subplot(2,2,1), c = fill(s, btw, 'g'); alpha(0.5);
    set(c,'EdgeColor','none')
    hold on
    plot(mean(diffMnFull),'-kx', 'MarkerSize', 8), set(gca,'xticklabel',{'1','2','4','8','16','32'}), title(['Mean, Full layer: ' chooseDatab ' ' chooseImCond]), xlabel('Density'), ylabel([texlabel(regressNm) 'symmetric - ' texlabel(regressNm) 'antisymm.'])
    line([1 6],[0 0],'LineStyle','--', 'Color',[0.5 0.5 0.5])
    axis square, xlim([1,6])

    ci = bootci(2000,@mean,diffMnRight);
    lower = ci(1,:); upper = ci(2,:);
    btw = [lower, fliplr(upper)];
    subplot(2,2,2), c = fill(s, btw, 'g'); alpha(0.5);
    set(c,'EdgeColor','none')
    hold on
    plot(mean(diffMnRight),'-kx', 'MarkerSize', 8), set(gca,'xticklabel',{'1','2','4','8','16','32'}), title(['Mean, Right half layer: ' chooseDatab ' ' chooseImCond]), xlabel('Density'), ylabel([texlabel(regressNm) 'symmetric - ' texlabel(regressNm) 'antisymm.'])
    line([1 6],[0 0],'LineStyle','--', 'Color',[0.5 0.5 0.5])
    axis square, xlim([1,6])

    ci = bootci(2000,@mean,diffVrFull);
    lower = ci(1,:); upper = ci(2,:);
    btw = [lower, fliplr(upper)];
    subplot(2,2,3), c = fill(s, btw, 'g'); alpha(0.5);
    set(c,'EdgeColor','none')
    hold on
    plot(mean(diffVrFull),'-kx', 'MarkerSize', 8), set(gca,'xticklabel',{'1','2','4','8','16','32'}), title(['Var, Full layer: ' chooseDatab ' ' chooseImCond]), xlabel('Density level'), ylabel([texlabel(regressNm) 'symmetric - ' texlabel(regressNm) 'antisymm.'])
    axis square, xlim([1,6])

    ci = bootci(2000,@mean,diffVrRight);
    lower = ci(1,:); upper = ci(2,:);
    btw = [lower, fliplr(upper)];
    subplot(2,2,4), c = fill(s, btw, 'g'); alpha(0.5);
    set(c,'EdgeColor','none')
    hold on
    plot(mean(diffVrRight),'-kx', 'MarkerSize', 8),  set(gca,'xticklabel',{'1','2','4','8','16','32'}), title(['Var, Right layer: ' chooseDatab ' ' chooseImCond]), xlabel('Density level'), ylabel([texlabel(regressNm) 'symmetric - ' texlabel(regressNm) 'antisymm.'])
    axis square, xlim([1,6])

    cd([BasisPath '/DATA/Figure7/noCM' filesep chooseGF filesep chooseDatab '_' chooseImCond filesep])
    if MODEL.regressFit == 1
        load('betadiffMean.mat')
        load('betadiffVariance.mat')
    elseif MODEL.regressFit == 2
        if exist('r2diffMean.mat') && exist('r2diffVariance.mat')
            load('r2diffMean.mat')
            load('r2diffVariance.mat')
        else
            error('NON EXISTING MAT-FILE(S) r2diffVariance.mat OR r2diffMean.mat not found')
        end
    end

    ci = bootci(2000,@mean,diffMnFull);
    lower = ci(1,:); upper = ci(2,:);
    btw = [lower, fliplr(upper)];
    subplot(2,2,1), c = fill(s, btw, 'y'); alpha(0.5);
    set(c,'EdgeColor','none')
    hold on
    plot(mean(diffMnFull),'-k')

    ci = bootci(2000,@mean,diffMnRight);
    lower = ci(1,:); upper = ci(2,:);
    btw = [lower, fliplr(upper)];
    subplot(2,2,2), c = fill(s, btw, 'y'); alpha(0.5);
    set(c,'EdgeColor','none')
    hold on
    plot(mean(diffMnRight),'-k')

    ci = bootci(2000,@mean,diffVrFull);
    lower = ci(1,:); upper = ci(2,:);
    btw = [lower, fliplr(upper)];
    subplot(2,2,3), c = fill(s, btw, 'y'); alpha(0.5);
    set(c,'EdgeColor','none')
    hold on
    plot(mean(diffVrFull),'-k')
    line([1 6],[0 0],'LineStyle','--', 'Color',[0.5 0.5 0.5])

    ci = bootci(2000,@mean,diffVrRight);
    lower = ci(1,:); upper = ci(2,:);
    btw = [lower, fliplr(upper)];
    subplot(2,2,4), c = fill(s, btw, 'y'); alpha(0.5);
    set(c,'EdgeColor','none')
    hold on
    plot(mean(diffVrRight),'-k')
    line([1 6],[0 0],'LineStyle','--', 'Color',[0.5 0.5 0.5])

%     print(fig, [Path2Save, 'Figure7c'], '-vector', '-depsc', '-r300')
    print(fig, [Path2Save, 'Figure7c'], '-depsc', '-r300')
    return
end

%% do RSA

if ~exist([BasisPath '/DATA/Figure8/'])
    mkdir([BasisPath '/DATA/Figure8'])
end

if ~exist([BasisPath '/DATA/Figure8/' chooseGF filesep chooseDatab '_', chooseImCond '.mat'])

    for id = 1:nFaces
        [corrCoef, corrDemCoef, eucCoef] = doRSA(PatternSaved, id, MODEL.nLayers, MODEL.density);

        corrCoefAll(:,:,:,id) = corrCoef;
        corrDemCoefAll(:,:,:,id) = corrDemCoef;
        eucCoefAll(:,:,:,id) = eucCoef;
    end


    if ~exist([BasisPath '/DATA/Figure8/' chooseGF])
        mkdir([BasisPath '/DATA/Figure8/' chooseGF])
    end
    cd ([BasisPath '/DATA/Figure8/' chooseGF])

    save([chooseDatab '_' chooseImCond '.mat'],'corrCoefAll','corrDemCoefAll','eucCoefAll')

else
    cd ([BasisPath '/DATA/Figure8/' chooseGF])
    load([chooseDatab '_' chooseImCond '.mat'],'corrCoefAll','corrDemCoefAll','eucCoefAll')
end

%% FIGURE 8: Plot RSA results, single density level

if MODEL.doFig8

    Xaxis = 1:MODEL.nLayers;
    s = [Xaxis, fliplr(Xaxis)];

    fig=figure('Position',[1 scrsz(4)/3 scrsz(3) scrsz(4)/3]);

    % corr
    subplot(1,3,1),
    mn = median(corrCoefAll,4);
    mn1 = mn(1,:,MODEL.pickDens); mn2 = mn(2,:,MODEL.pickDens);
    upper= prctile(corrCoefAll,75,4);
    lowper= prctile(corrCoefAll,25,4);
    top1 = upper(1,:,MODEL.pickDens); bot1 = lowper(1,:,MODEL.pickDens);
    btw1 = [top1, fliplr(bot1)];
    c = fill(s, btw1, 'r'); alpha(0.2);
    set(c,'EdgeColor','none')
    hold on;
    top2 = upper(2,:,MODEL.pickDens); bot2 = lowper(2,:,MODEL.pickDens);
    btw2 = [top2, fliplr(bot2)];
    d = fill(s, btw2, 'b'); alpha(0.2);
    set(d,'EdgeColor','none')
    plot(Xaxis, mn1, 'r', Xaxis, mn2, 'b', 'LineWidth', 2), xlabel('Layer'), ylabel('Correlation (Spearman)'), title(['Corr, ' chooseDatab ' ' chooseImCond ', dens: ' num2str(MODEL.density(MODEL.pickDens))])

    % corr right dem
    subplot(1,3,2),
    mn = median(corrDemCoefAll,4);
    mn1 = mn(1,:,MODEL.pickDens); mn2 = mn(2,:,MODEL.pickDens);
    upper= prctile(corrDemCoefAll,75,4);
    lowper= prctile(corrDemCoefAll,25,4);
    top1 = upper(1,:,MODEL.pickDens); bot1 = lowper(1,:,MODEL.pickDens);
    btw1 = [top1, fliplr(bot1)];
    c = fill(s, btw1, 'r'); alpha(0.2);
    set(c,'EdgeColor','none')
    hold on;
    top2 = upper(2,:,MODEL.pickDens); bot2 = lowper(2,:,MODEL.pickDens);
    btw2 = [top2, fliplr(bot2)];
    d = fill(s, btw2, 'b'); alpha(0.2);
    set(d,'EdgeColor','none')
    plot(Xaxis, mn1, 'r', Xaxis, mn2, 'b', 'LineWidth', 2), xlabel('Layer'), ylabel('Correlation (Spearman)'), title(['Corr Dem, ' chooseDatab ' ' chooseImCond ', dens: ' num2str(MODEL.density(MODEL.pickDens))])

    % euc right
    subplot(1,3,3),
    mn = median(eucCoefAll,4);
    mn1 = mn(1,:,MODEL.pickDens); mn2 = mn(2,:,MODEL.pickDens);
    upper= prctile(eucCoefAll,75,4);
    lowper= prctile(eucCoefAll,25,4);
    top1 = upper(1,:,MODEL.pickDens); bot1 = lowper(1,:,MODEL.pickDens);
    btw1 = [top1, fliplr(bot1)];
    c = fill(s, btw1, 'r'); alpha(0.2);
    set(c,'EdgeColor','none')
    hold on;
    top2 = upper(2,:,MODEL.pickDens); bot2 = lowper(2,:,MODEL.pickDens);
    btw2 = [top2, fliplr(bot2)];
    d = fill(s, btw2, 'b'); alpha(0.2);
    set(d,'EdgeColor','none')
    plot(Xaxis, mn1, 'r', Xaxis, mn2, 'b', 'LineWidth', 2), xlabel('Layer'), ylabel('Correlation (Spearman)'), title(['Euc, ' chooseDatab ' ' chooseImCond ', dens: ' num2str(MODEL.density(MODEL.pickDens))])

%     print(fig,[Path2Save, 'Figure8'], '-vector', '-dsvg', '-r300');
    print(fig,[Path2Save, 'Figure8'], '-dsvg', '-r300');

    return
end

%% FIGURE 9: Plot t-test of RSA results, all density levels

if MODEL.doFig9

    fig=figure('Position',[1 scrsz(4)/3 scrsz(3) scrsz(4)/3]);

    % corr
    for j = 1:length(MODEL.density)
        for i = 1:MODEL.nLayers
            [h,p,ci,stats] = ttest((corrCoefAll(1,i,j,:)),(corrCoefAll(2,i,j,:))); %#ok<*ASGLU> %1 is viewpoint model, 2 is symmetry model: ttest compares corr values between both models across all IDs/blocks
            tvalCorr(j,i) = stats.tstat;
        end
    end

    subplot(1,3,1),
    s = surf(tvalCorr); view(0,90);
    xlabel('Layer'), ylabel('Density'),
    title(['Corr, ' chooseDatab ' ' chooseImCond]),
    colorbar, colormap(redblue)
    yticks([1 2 3 4 5 6])
    yticklabels({'1','2','4','8','16','32'})
    set(s,'linestyle','none')
    s.FaceColor = 'interp';
    maxT = max(tvalCorr(:)); minT = min(tvalCorr(:));
    if abs(maxT) > abs(minT)
        caxis([-(abs(maxT)) abs(maxT)]);
    elseif abs(minT) > abs(maxT)
        caxis([-(abs(minT)) abs(minT)]);
    else
        caxis([-(abs(minT)) abs(minT)]);
    end

    for j = 1:length(MODEL.density)
        for i = 1:MODEL.nLayers
            [h,p,ci,stats] = ttest((corrDemCoefAll(1,i,j,:)),(corrDemCoefAll(2,i,j,:)));
            tvalCorrDem(j,i) = stats.tstat;
        end
    end

    subplot(1,3,2),
    s = surf(tvalCorrDem); view(0,90);
    xlabel('Layer'), ylabel('Density')
    title(['Corr Dem, ' chooseDatab ' ' chooseImCond]),
    colorbar, colormap(redblue)
    yticks([1 2 3 4 5 6])
    yticklabels({'1','2','4','8','16','32'})
    set(s,'linestyle','none')
    s.FaceColor = 'interp';
    maxT = max(tvalCorrDem(:)); minT = min(tvalCorrDem(:));
    if abs(maxT) > abs(minT)
        caxis([-(abs(maxT)) abs(maxT)]); %#ok<*CAXIS>
    elseif abs(minT) > abs(maxT)
        caxis([-(abs(minT)) abs(minT)]);
    else
        caxis([-(abs(minT)) abs(minT)]);
    end

    for j = 1:length(MODEL.density)
        for i = 1:MODEL.nLayers
            [h,p,ci,stats] = ttest((eucCoefAll(1,i,j,:)),(eucCoefAll(2,i,j,:)));
            tvalEuc(j,i) = stats.tstat;
        end
    end

    subplot(1,3,3),
    s = surf(tvalEuc); view(0,90);
    xlabel('Layer'), ylabel('Density')
    title(['Euc, ' chooseDatab ' ' chooseImCond]),
    colorbar, colormap(redblue)
    yticks([1 2 3 4 5 6])
    yticklabels({'1','2','4','8','16','32'})
    set(s,'linestyle','none')
    s.FaceColor = 'interp';
    maxT = max(tvalEuc(:)); minT = min(tvalEuc(:));
    if abs(maxT) > abs(minT)
        caxis([-(abs(maxT)) abs(maxT)]);
    elseif abs(minT) > abs(maxT)
        caxis([-(abs(minT)) abs(minT)]);
    else
        caxis([-(abs(minT)) abs(minT)]);
    end

%     print(fig,[Path2Save, 'Figure9'], '-vector', '-depsc', '-r300');
    print(fig,[Path2Save, 'Figure9'], '-depsc', '-r300');

    return
end


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
% checkSumR2 = sum(R2All);

end





