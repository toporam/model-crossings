%% MAIN: model-crossings
clear, clc 

%-------------------------------------------------------------------------
% GOAL: replicate simulations and show figures
%-------------------------------------------------------------------------
% Set path

%BasisPath = % <-- Set this to directory containing ImageDatasets folder
%path2code = % <-- Set this to directory containing all code

BasisPath = '/Users/ramirezfm/Documents/PROJECTS/crossings/model-crossings-dir-1';
path2code = '/Users/ramirezfm/Documents/GitHub/model-crossings';

% ----------------------------------------
% Before running this script, you need to:
% ----------------------------------------
% 1. Obtain permission and download the KDEF and Radboud face databases
% 2. Create directory BasisPath (path specified in line 14)
% 3. Create a sub-directory at BasisPath named rawImg
% 4. Place in 'BasisPath/rawImg' the two directories downloaded in step 1, 
% above (one containing the KDEF database, the other the Radboud database)
% 5. Make sure that the paths specified in lines 33 and 34, below, match the
% location of the KDEF and Radboud databases
% ----------------------------------------

cd('/')
addpath(genpath(path2code),'-end')

kdefPath = [BasisPath '/rawImg/KDEF_and_AKDEF/KDEF']; %<-- set path where kdef directory is downloaded 
rafdPath = [BasisPath '/rawImg/RafDDownload-F848167628A5D2F935A500BB316AD773']; %<-- set path where rafd directory is downloaded 

% prep downloaded image databases for use in the model:
prepDatabImages(BasisPath, kdefPath, rafdPath) 

cd(BasisPath)
if ~exist([BasisPath '/Results']) %#ok<*EXIST>
    mkdir('Results')
end

if ~exist([BasisPath '/DATA'])
    mkdir('DATA') 
end
%-------------------------------------------------------------------------
% MODEL TYPES: a structure for each of the 4 model variants is saved with 
% its defining parameters
%
% MODEL1: radboud/pixel
% MODEL2: radboud/s1
% MODEL3: kdef/pixel
% MODEL4: kdef/s1

%-------------------------------------------------------------------------
% SET MODEL1 PARAMETERS AND FLAGS
%-------------------------------------------------------------------------
% Specify parameters of the Feedforward Randomly-Connected 2-Hemishpere
% (RC-2H) networks to be generated
%-------------------------------------------------------------------------
MODEL1.mkConnect = 1; % [0/1; no/yes] 1 will generate network connections (necessary for first time running the model). 0 will load pre-generated network if it exists
    MODEL1.nLayers = 8;
    MODEL1.nUnitsL1 = 2^11; MODEL1.nUnitsRest = 1024; % nUnitsL1= # units for HALF layer 1, nUnitsRest= # units for each full layer 2-8
    MODEL1.density = [1,2,4,8,16,32];
    MODEL1.biProb = 0.98; % Binomial probability of a projection crossing into the ispilateral hemisphere for projections from Layer 1 to Layer 2
    MODEL1.k = 4; % Logistic growth rate

%-------------------------------------------------------------------------
MODEL1.applyGF = 1; %[0/1] 1 to apply gain field to activation patterns. 0 to have no gain field
%-------------------------------------------------------------------------
% Show face images as input to the network
%-------------------------------------------------------------------------
MODEL1.seeImg = 0; % [0/1] 1 will display each face identity at all 5 orientations-- warning: will produce lots of figures!
%-------------------------------------------------------------------------
%**************************************************************************
% Define image database to be used as input to the generated networks, as
% well as the model variant {pixel/s1}
%**************************************************************************

% Image database
% **************
MODEL1.indDatab = 1; % 1= radboud, 2= kdef

% Model variant
% **************
MODEL1.indImCond = 1; %1= pixel, 2= s1

% Do Cortical Magnification
% **************
MODEL1.indCortMag = 1; %[1/2] 1=include CM (default), 2=don't include CM

%-------------------------------------------------------------------------
% SET MODEL2 PARAMETERS AND FLAGS
%-------------------------------------------------------------------------
MODEL2 = MODEL1;
MODEL2.indImCond = 2; % 1= pixel, 2= s1

MODEL2.seeMadeImg = 0; % [0/1] if making S1 representations, 1 = view images as they are being produced -- warning: will produce lots of figures!

%-------------------------------------------------------------------------
% SET MODEL3 PARAMETERS AND FLAGS
%-------------------------------------------------------------------------
MODEL3 = MODEL1;
MODEL3.indDatab = 2; %1= radboud, 2= kdef

%-------------------------------------------------------------------------
% SET MODEL4 PARAMETERS AND FLAGS
%-------------------------------------------------------------------------
MODEL4 = MODEL1;
MODEL4.indImCond = 2; % %1= pixel, 2= s1
MODEL4.indDatab = 2; % 1= radboud, 2= kdef

MODEL4.seeMadeImg = 0; 

%% Make FIGURE 4 
MODEL1.doFig4 = 1; % [0/1] display the intermediate and final probability distributions. By default, will save Figure 3.

callCrossingsModel(MODEL1, BasisPath)

%% Make FIGURE 5
MODEL1.barColor = 0.9;
MODEL1.doFig6 = 0;
MODEL1.doFig5 = 1; % [0/1] plot median of distributions of mean and variance of pixel intensities for whole images. By default, will save Figure 4.
MODEL1.doFig4 = 0; 
callCrossingsModel(MODEL1, BasisPath) 

MODEL2.barColor = 0.1;
MODEL2.doFig6 = 0;
MODEL2.doFig5 = 1;
MODEL2.doFig4 = 0;
callCrossingsModel(MODEL2, BasisPath) 

MODEL3.barColor = 0.9;
MODEL3.doFig6 = 0;
MODEL3.doFig5 = 1;
MODEL3.doFig4 = 0; 
callCrossingsModel(MODEL3, BasisPath) 

MODEL4.barColor = 0.1;
MODEL4.doFig6 = 0;
MODEL4.doFig5 = 1;
MODEL4.doFig4 = 0; 
callCrossingsModel(MODEL4, BasisPath) 

%% Make FIGURE 6
MODEL1.barColor = 0.9;
MODEL1.doFig6 = 1; % [0/1] plot median of distributions of mean and variance of pixel intensities for half images. By default, will save Figure 5.
MODEL1.doFig5 = 0; 
MODEL1.doFig4 = 0;
callCrossingsModel(MODEL1, BasisPath) 

MODEL2.barColor = 0.1;
MODEL2.doFig6 = 1;
MODEL2.doFig5 = 0; 
MODEL2.doFig4 = 0;
callCrossingsModel(MODEL2, BasisPath) 

MODEL3.barColor = 0.9;
MODEL3.doFig6 = 1;
MODEL3.doFig5 = 0; 
MODEL3.doFig4 = 0;
callCrossingsModel(MODEL3, BasisPath) 

MODEL4.barColor = 0.1;
MODEL4.doFig6 = 1;
MODEL4.doFig5 = 0; 
MODEL4.doFig4 = 0;
callCrossingsModel(MODEL4, BasisPath) 

%% Make FIGURE 7

MODEL1.indCortMag = 1; %[1/2] 1=include CM, 2=don't include CM
MODEL1.pickLay = 1; % [1:8] layer # for activation pattern analysis
MODEL1.pickDens = 1; % [1:6 = density of 2^0:2^5] density-level for computing the mean and variance of activation patterns
MODEL1.regressFit = 2; % [1/2] which measure to use to describe fit of variance/norm trends with regressors; 1 is beta diffs, 2 is R^2 diffs 
MODEL1.doFig7 = 1; % [0/1] plot median of distributions of variance of activations for layer 1, and compare CM and no-CM condition. By default, will save Figure 6.
MODEL1.doFig6 = 0; MODEL1.doFig5 = 0; MODEL1.doFig4 = 0;
callCrossingsModel(MODEL1, BasisPath) 

MODEL5 = MODEL1;
MODEL5.indCortMag = 2;  % same as model 1 except don't include CM 
callCrossingsModel(MODEL5,  BasisPath) 

MODEL2.indCortMag = 1;
MODEL2.pickLay = 1; 
MODEL2.pickDens = 1;
MODEL2.regressFit = 2;
MODEL2.doFig7 = 1;
MODEL2.doFig6 = 0; MODEL2.doFig5 = 0; MODEL2.doFig4 = 0;   
callCrossingsModel(MODEL2,  BasisPath) 

MODEL6 = MODEL2; 
MODEL6.indCortMag = 2; 
callCrossingsModel(MODEL6,  BasisPath) 

MODEL3.indCortMag = 1;
MODEL3.pickLay = 1; 
MODEL3.pickDens = 1;
MODEL3.regressFit = 2;
MODEL3.doFig7 = 1; 
MODEL3.doFig6 = 0; MODEL3.doFig5 = 0; MODEL3.doFig4 = 0; 
callCrossingsModel(MODEL3,  BasisPath) 

MODEL7 = MODEL3;
MODEL7.indCortMag = 2; 
callCrossingsModel(MODEL7,  BasisPath) 

MODEL4.indCortMag = 1;
MODEL4.pickLay = 1; 
MODEL4.pickDens = 1;
MODEL4.regressFit = 2;
MODEL4.doFig7 = 1; 
MODEL4.doFig6 = 0; MODEL4.doFig5 = 0; MODEL4.doFig4 = 0; 
callCrossingsModel(MODEL4, BasisPath)

MODEL8 = MODEL4;
MODEL8.indCortMag = 2; 
callCrossingsModel(MODEL8, BasisPath) 

%% Make FIGURE 8
MODEL1.pickDens = 5;
MODEL1.doFig8 = 1; % [0/1] plot RSA results for a single density level. By default, will save Figure 7. 
MODEL1.doFig7 = 0; MODEL1.doFig6 = 0; MODEL1.doFig5 = 0; MODEL1.doFig4 = 0; 
callCrossingsModel(MODEL1, BasisPath) 

MODEL2.pickDens = 5;
MODEL2.doFig8 = 1;
MODEL2.doFig7 = 0; MODEL2.doFig6 = 0; MODEL2.doFig5 = 0; MODEL2.doFig4 = 0; 
callCrossingsModel(MODEL2, BasisPath) 

MODEL3.pickDens = 5;
MODEL3.doFig8 = 1; 
MODEL3.doFig7 = 0; MODEL3.doFig6 = 0; MODEL3.doFig5 = 0; MODEL3.doFig4 = 0; 
callCrossingsModel(MODEL3, BasisPath) 

MODEL4.pickDens = 5;
MODEL4.doFig8 = 1; 
MODEL4.doFig7 = 0; MODEL4.doFig6 = 0; MODEL4.doFig5 = 0; MODEL4.doFig4 = 0; 
callCrossingsModel(MODEL4, BasisPath) 

%% Make FIGURE 9
MODEL1.doFig9 = 1; % [0/1] plot RSA results for all density levels. By default, will save Figure 8. 
MODEL1.doFig8 = 0; MODEL1.doFig7 = 0; MODEL1.doFig6 = 0; MODEL1.doFig5 = 0; MODEL1.doFig4 = 0; 
callCrossingsModel(MODEL1, BasisPath) 
auxTime = clock; display(['Time: ' num2str(auxTime(4)), ':', num2str(auxTime(5)) ' hrs']),% beep

MODEL2.doFig9 = 1;
MODEL2.doFig8 = 0; MODEL2.doFig7 = 0; MODEL2.doFig6 = 0; MODEL2.doFig5 = 0; MODEL2.doFig4 = 0; 
callCrossingsModel(MODEL2, BasisPath) 
auxTime = clock; display(['Time: ' num2str(auxTime(4)), ':', num2str(auxTime(5)) ' hrs']),% beep

MODEL3.doFig9 = 1; 
MODEL3.doFig8 = 0; MODEL3.doFig7 = 0; MODEL3.doFig6 = 0; MODEL3.doFig5 = 0; MODEL3.doFig4 = 0; 
callCrossingsModel(MODEL3, BasisPath) 
auxTime = clock; display(['Time: ' num2str(auxTime(4)), ':', num2str(auxTime(5)) ' hrs']), % beep

MODEL4.doFig9 = 1; 
MODEL4.doFig8 = 0; MODEL4.doFig7 = 0; MODEL4.doFig6 = 0; MODEL4.doFig5 = 0; MODEL4.doFig4 = 0; 
callCrossingsModel(MODEL4, BasisPath) 
auxTime = clock; display(['Time: ' num2str(auxTime(4)), ':', num2str(auxTime(5)) ' hrs']), % beep

%% 