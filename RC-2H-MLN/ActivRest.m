function RestActisC = ActivRest(indxLRest,L1Actis,nLayers,density)

%layer 2
lay=1; 
p=0;
for d = density
    p = p+1; 
    for indDens = 1:d   
        temp(:,:,indDens) = L1Actis(indxLRest(:,indDens,lay),:,p); %gets the activation values of layer 1 untis according to the indices of indxLRest   
    end
    RestActis(:,:,p) = mean(temp,3); %averages to get layer 2 values
end

RestActisC{lay} = RestActis;


%layers 3-8

for lay = 2:(nLayers-1)
    p=0;
    for d = density
        p = p+1;
        for indDens = 1:d   
            aux(:,:,:) = RestActisC{lay-1};
            temp(:,:,indDens) = aux(indxLRest(:,indDens,lay),:,p); %gets the activation values of the previous layer according to the indices of the current layer        
        end
        RestActis(:,:,p) = mean(temp,3);
    end
    
RestActisC{lay} = RestActis; %RestActisC is a cell of size 7, for layers 2-8
end 

