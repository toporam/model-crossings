function S1Img = getS1(chooseDatab, nFaces, BasisPath, seeMadeImg, imgAsIs)

% December 2014, adapted July 2021

% Obtain S1 representation of each image based on HMAX model as implemented
% by Serre et al., 2004.

patchSizes = [4 8 12 16]; 
			  
numPatchSizes = length(patchSizes);

view = [0,45,90,135,180];

s=0;
for id = 1:nFaces
    for ori = 1:length(view)
        s = s+1;
        if imgAsIs == 0 
            aux = imread([BasisPath '/ImageDatasets/pixel/' chooseDatab filesep chooseDatab num2str(view(ori)) '_' num2str(id) '.bmp']);
        elseif imgAsIs == 1    
            aux = imread([BasisPath '/ImageDatasets/pixel/' chooseDatab 'Only' filesep chooseDatab num2str(view(ori)) '_' num2str(id) '.bmp']);
        end
        cI{s,1} = double(aux);
    end
end

if isempty(cI{1}) 
  error('No training images were loaded -- did you remember to change the path names?');
end
  
% Below the c1 prototypes are extracted from the images read from file.
% Although not used here, needed as input for function that extracts s1.

fprintf('reading patches');
cPatches = load('PatchesFromNaturalImages250per4sizes','cPatches');
cPatches = cPatches.cPatches;


%---------- defauls ---------
rot = [90 -45 0 45];
c1ScaleSS = [1:2:18]; %#ok<NBRAK>
RF_siz    = [7:2:39]; %#ok<NBRAK> 
c1SpaceSS = [8:2:22]; %#ok<NBRAK>
minFS     = 7;
maxFS     = 39;
div = [4:-.05:3.2]; %#ok<NBRAK>
Div       = div;
%-----------------------------

fprintf(1,'Initializing gabor filters -- full set...');
% Creates the gabor filters use to extract the S1 layer
[fSiz,filters,c1OL] = init_gabor_cos(rot, RF_siz, Div);

fprintf(1,'done\n');

% The actual S1 features are computed below for each one of the training/testing directories
tic
S1res = extractS1forcell(filters,fSiz,c1SpaceSS,c1ScaleSS,c1OL,cPatches,cI,numPatchSizes); %#ok<SAGROW>
toc


S1Img = S1addOutput(cI, S1res, nFaces, seeMadeImg);

