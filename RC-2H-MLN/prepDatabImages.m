function prepDatabImages(BasisPath, kdefPath, rafdPath)

% This function formats images from the Radboud Faces Database and
% Karolinska Directed Emotional Faces (KDEF) dataset as required for their
% use in the present model. Before running function, you must have
% downloaded both databases and located them in the relevant directories as
% specified in mainMkFigs.m.

%% Make radboudOnly directory
    
pathradboudOnly = [BasisPath '/ImageDatasets/pixel/radboudOnly'];

if ~exist(pathradboudOnly,'dir')
    mkdir(pathradboudOnly)

    views = {'000', '045', '090', '135', '180'};
    views2 = [0, 45, 90, 135, 180];
    
    cd(rafdPath)
    
    for j = 1:5
        S = dir(fullfile(rafdPath, ['*', views{j}, '*ca*neutral_frontal.jpg'])); % get all neutral, frontal gaze adult faces
        N = setdiff({S.name},{'.','..'});
    
        if numel(N) ~= 57
            warning('Number of IDs is incorrect-- should have 57 IDs')
        end
    
        p = 0;
        for i = 1:numel(N)
            p = p + 1;
            img = imread(N{i});
            grImg = rgb2gray(img);
            grImg = imresize(grImg,[594 395]);
            grImg = grImg(25:465, :); %resize and crop so that longer side is 441 pixels
            imwrite(grImg, [pathradboudOnly filesep 'radboud' num2str(views2(j)) '_' num2str(p) '.bmp'])
        end 
    end
end 

%% Make radboud directory 

% Resize images to ~ 8.7° visual angle and center on gray square background 

pathradboud = [BasisPath '/ImageDatasets/pixel/radboud'];

if ~exist(pathradboud,'dir')
    mkdir(pathradboud)

    S = dir(fullfile(pathradboudOnly,'*'));
    N = setdiff({S.name},{'.','..'});
    
    
    for i = 1:numel(N)
        cd(pathradboudOnly)
        graySq = ones(441,441)*128;    
        img = imread(N{i});
        img = imresize(img,[369 331]);
        graySq(37:405, 56:386) = img;
        imwrite(uint8(graySq), [pathradboud filesep N{i}])
    end 
end 

%% Make kdefOnly directory
% select neutral face images from the KDEF database, using the "B" version
% as default. It converts images to grayscale/ resizes, and saves in new
% directory

pathkdefOnly = [BasisPath '/ImageDatasets/pixel/kdefOnly'];

if ~exist(pathkdefOnly,'dir')
    mkdir(pathkdefOnly)

    cd(kdefPath)
    addpath(genpath(kdefPath),'-end')

    S = dir(fullfile(kdefPath,'B*'));
    N = setdiff({S([S.isdir]).name},{'.','..'});
    
    for i = 1:numel(N)
        img = imread([N{i}, 'NEFR.JPG']);
        grImg = rgb2gray(img);
        lumB(i,1) = mean(grImg(:)); %#ok<*AGROW>
        img = imread(['A' N{i}(2:end), 'NEFR.JPG']);
        grImg = rgb2gray(img);
        lumA(i,1) = mean(grImg(:));
    
        img = imread([N{i}, 'NEHR.JPG']);
        grImg = rgb2gray(img);    
        lumB(i,2) = mean(grImg(:));
        img = imread(['A' N{i}(2:end), 'NEHR.JPG']);
        grImg = rgb2gray(img);
        lumA(i,2) = mean(grImg(:));
    
        img = imread([N{i}, 'NES.JPG']);
        grImg = rgb2gray(img);    
        lumB(i,3) = mean(grImg(:));
        img = imread(['A' N{i}(2:end), 'NES.JPG']);
        grImg = rgb2gray(img);
        lumA(i,3) = mean(grImg(:));
    
        img = imread([N{i}, 'NEHL.JPG']);
        grImg = rgb2gray(img);    
        lumB(i,4) = mean(grImg(:));
        img = imread(['A' N{i}(2:end), 'NEHL.JPG']);
        grImg = rgb2gray(img);
        lumA(i,4) = mean(grImg(:));
    
        img = imread([N{i}, 'NEFL.JPG']);
        grImg = rgb2gray(img);    
        lumB(i,5) = mean(grImg(:));
        img = imread(['A' N{i}(2:end), 'NEFL.JPG']);
        grImg = rgb2gray(img);
        lumA(i,5) = mean(grImg(:));
    end 
    
    p = 0;
    for i = 1:numel(N)
    
        % IDs containing images where mean luminance is more than 4 SD away 
        % from mean of database are replaced with images of same ID from set A. 
        % If neither version is usable, ID is excluded
    
        if any(lumB(i,:) > (mean(lumB(:)) + 4*std(lumB(:)))) ||  any(lumB(i,:) < (mean(lumB(:)) - 4*std(lumB(:))))
            if any(lumA(i,:) > (mean(lumA(:)) + 4*std(lumA(:)))) ||  any(lumA(i,:) < (mean(lumA(:)) - 4*std(lumA(:))))
                continue
            else 
                st = ['A' N{i}(2:end)];
            end
        else 
             st = N{i};
        end 
    
        p = p+1;
    
        img = imread([st, 'NEFR.JPG']);
        grImg = rgb2gray(img);
        grImg = imresize(grImg,[441 326]); %resize so that longer side is 441 pixels
        imwrite(grImg, [pathkdefOnly filesep 'kdef0_' num2str(p) '.bmp'])
    
        img = imread([st, 'NEHR.JPG']);
        grImg = rgb2gray(img);
        grImg = imresize(grImg,[441 326]);
        imwrite(grImg, [pathkdefOnly filesep 'kdef45_' num2str(p) '.bmp'])
    
        img = imread([st, 'NES.JPG']);
        grImg = rgb2gray(img);
        grImg = imresize(grImg,[441 326]);
        imwrite(grImg, [pathkdefOnly filesep 'kdef90_' num2str(p) '.bmp'])
    
        img = imread([st, 'NEHL.JPG']);
        grImg = rgb2gray(img);
        grImg = imresize(grImg,[441 326]);
        imwrite(grImg, [pathkdefOnly filesep 'kdef135_' num2str(p) '.bmp'])
    
        img = imread([st, 'NEFL.JPG']);
        grImg = rgb2gray(img);
        grImg = imresize(grImg,[441 326]);
        imwrite(grImg, [pathkdefOnly filesep 'kdef180_' num2str(p) '.bmp'])
    
    end
    
    if p ~= 68
        warning('Number of IDs is incorrect-- should have 68 IDs')
    end 
    
    clear img grImg
end 

%% Make kdef directory
% Resize images to ~ 8.7° visual angle and center on gray square background 

pathkdef = [BasisPath '/ImageDatasets/pixel/kdef'];

if ~exist(pathkdef,'dir')
    mkdir(pathkdef)
    
    S = dir(fullfile(pathkdefOnly,'*'));
    N = setdiff({S.name},{'.','..','.DS_Store'});
    
    
    for i = 1:numel(N)
        cd(pathkdefOnly)
        graySq = ones(441,441)*128;    
        img = imread(N{i});
        img = imresize(img,[341 253]);
        graySq(51:391, 95:347) = img;
        imwrite(uint8(graySq), [pathkdef filesep N{i}])
    end 
    
    clear img graySq 

end 