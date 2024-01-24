function [L1Actis, RestActisC] = defineActiv(imageV, indxL1, indxLRest, nLayers, density)
          
% Layer 1
p=0;
for d = density
    p = p+1;
    for indDens = 1:d   
        temp(:,:,indDens) = imageV(indxL1(:,indDens),:); % gets the values at image locations according to indxL1, and following density value
    end 
    L1Actis(:,:,p) = mean(temp,3); %if density is higher than 1, average the image values
end 


% Layer 2
RestActisC = ActivRest(indxLRest,L1Actis,nLayers,density);
