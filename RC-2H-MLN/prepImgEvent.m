function imageV = prepImgEvent(chooseDatab, chooseImCond, id, seeImg, BasisPath)
% For the unique combination of database (KDEF, RaFD), image representation
% (pixel, S1), and face identity (id) indicated by the input parameters,
% function will read the images associated with the five face views of that
% id and provide them as output in matrix form, ordered as defined by the
% vector views, defined below.
%
% chooseDatab:  choose database ('pixel' or 's1')
%
% chooseImCond: choose image representation ('pixel' or 's1')
%
% id: index of specific face ID to be loaded
%
% seeImage: flag indicating whether to show the 5 views of that id or not
%
% BasisPath: path2dir with sub-directories for both databases
%
% ImageV: matrix with dimensions (441*441) by 5--i.e., nPix^2*nView
%
% Cambria Revsine and Fernando Ramirez, Feb 2023


views = [0, 45, 90, 135, 180];

if strcmp(chooseImCond, 'pixel')    
    for ori = 1:length(views)
        imageS = imread([BasisPath '/ImageDatasets/pixel/' chooseDatab filesep chooseDatab num2str(views(ori)) '_' num2str(id) '.bmp']);
        imageV(:,ori) = imageS(:);
    end
elseif strcmp(chooseImCond, 's1')
    load([BasisPath '/ImageDatasets/s1/s1Img_' chooseDatab '.mat'], 'S1Img');
    for ori = 1:length(views)
        imageS = cell2mat(S1Img(id,ori)); 
        imageV(:,ori) = imageS(:);
    end        
end

if seeImg
    scrsz = get(groot, 'Screensize');
    fig = figure('Position', [1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
    for ori = 1:length(views)
        imageS = reshape(imageV(:,ori),[441,441]);
        subplot(1,length(chooseViews),ori), imagesc(imageS), colormap gray, axis square, set(gca,'xtick',[], 'ytick',[])
    end
end

