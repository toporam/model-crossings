function cumS1 = extractS1forcell(filters,fSiz,c1SpaceSS,c1ScaleSS,c1OL,cPatches,cImages,numPatchSizes)
%this function is a wrapper of C1. For each image in the cell cImages, 
%it extracts all the values of the C1 and S1 layers for all the prototypes
%in the cell cPatches. The result cumS1 is a matrix of size ?

%total_number_of_patches \times number_of_images where
%total_number_of_patches is the sum over i = 1:numPatchSizes of
%length(cPatches{i}) and number_of_images is length(cImages) The C1
%parameters used are given as the variables
%filters,fSiz,c1SpaceSS,c1ScaleSS,c1OL for more detail regarding these
%parameters see the help entry for C1
%
%See also C1

%% a bug was fixed on Jul 01 2005

numPatchSizes = min(numPatchSizes,length(cPatches));
%all the patches are being flipped. This is becuase in matlab conv2 is much faster than filter2
for i = 1:numPatchSizes
  [siz,numpatch] = size(cPatches{i});
  siz = sqrt(siz/4);
  for j = 1:numpatch
    tmp = reshape(cPatches{i}(:,j),[siz,siz,4]);
    tmp = tmp(end:-1:1,end:-1:1,:);
    cPatches{i}(:,j) = tmp(:);
  end
end

% fdo replace loop c2 for c1 & save c1 & s1

cumC1 = [];
cumS1 = [];
  
for i = 1:length(cImages) %for every input image
    
  fprintf(1,'%d:',i);
  stim = cImages{i};
  img_siz = size(stim);

  [c1,s1] = C1(stim, filters, fSiz, c1SpaceSS, c1ScaleSS, c1OL);
  
  cumC1 = [cumC1; c1];
  cumS1 = [cumS1; s1];
    
end

