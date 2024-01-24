function [corrCoef, corrDemCoef, eucCoef] = doRSA(PatternSaved, id, nLayers, density)

% ALL RESULTS FOR RIGHT OF NETWORK = LEFT OF IMAGE (& vice versa) -----------

scrsz = get(groot, 'Screensize');

% Compute Corr and Euc distance (and demean)
for i = 1:nLayers
    for j = 1:length(density)
        pat = PatternSaved{id}{1,i};
        temp = pat(:,:,j);
        auxUnits = length(temp)/2;
        tempRight = temp(auxUnits+1:end,:); % 1 patterns vector for each orientation, from right hemisphere
        tempRightDem = tempRight-repmat(mean(tempRight,2),1,5); %data demean: subtract each unit's average response to all conditions from each unit at each condition
     
        corrRight{i}(:,:,j) = squareform(pdist(tempRight','Correlation'));
        corrRightV(:,j,i) = pdist(tempRight','Correlation'); % vector form

        corrRightDem{i}(:,:,j) = squareform(pdist(tempRightDem','Correlation'));
        corrRightDemV(:,j,i) = pdist(tempRightDem','Correlation');

        eucRight{i}(:,:,j) = squareform(pdist(tempRight','Euclidean'));
        eucRightV(:,j,i) = pdist(tempRight','Euclidean');
    end
end
    
clear pat
clear temp
clear tempRight
clear tempRightDem

%% RSA

% define model templates: only upper triangle of models

Viewpoint = [0.25 0.5 0.75 1 0.25 0.5 0.75 0.25 0.5 0.25];
Symmetry =[0.5 1 0.5 0 0.5 0 0.5 0.5 1 0.5];

% corr
for j = 1:length(density)
    for i = 1:nLayers
        corrCoef(1,i,j) = corr(Viewpoint',corrRightV(:,j,i),'type','Spearman'); 
        corrCoef(2,i,j) = corr(Symmetry',corrRightV(:,j,i),'type','Spearman');
    end
end

% corr dem 
for j = 1:length(density)
    for i = 1:nLayers
        corrDemCoef(1,i,j) = corr(Viewpoint',corrRightDemV(:,j,i),'type','Spearman'); 
        corrDemCoef(2,i,j) = corr(Symmetry',corrRightDemV(:,j,i),'type','Spearman');
    end
end

% euc
for j = 1:length(density)
    for i = 1:nLayers
        eucCoef(1,i,j) = corr(Viewpoint',eucRightV(:,j,i),'type','Spearman'); 
        eucCoef(2,i,j) = corr(Symmetry',eucRightV(:,j,i),'type','Spearman'); 
    end
end
