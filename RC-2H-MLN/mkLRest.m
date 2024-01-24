function indxLRest = mkLRest(nUnitsL1, nUnitsRest, density,nLayers, biProb)

distStep = (biProb-0.5)/(nLayers-2); % value the binomial probability decreases by in each layer

indxLRest = NaN(nUnitsRest,max(density),nLayers-1);

% layer 2
lay=1;

for i = 1:nUnitsRest
    for d = 1:max(density)
        if i <= nUnitsRest/2
            biVal(d) = binornd(1,1-biProb); % selects from binomial distribution to get an output of either 0 or 1, d times. biProb starting high (0.98) will almost always give a 0 for the first half of layer 2...
        elseif i > nUnitsRest/2
            biVal(d) = binornd(1,biProb); % but will almost always give a 1 for the 2nd half of L2
        end
    end
    biValInt = biVal*nUnitsL1; % multiply vector of 0s and 1s by half the length of L1 (nUnitsL1 = half the length)
    indW = randperm(nUnitsL1,max(density)); %get d random indices from half length of L1
    indxLRest(i,:,lay) = indW+biValInt; %add to biValInt to get final indices; when biValInt is 0, index will come from left half of L1, when biValInt = nUnitsL1, index will come from right half of L1
    clear biVal biValInt indW
end

% layers 3-8

for lay = 2:(nLayers-1)
    biProb = biProb - distStep; % decrease biProb for every layer- ends at 0.5
    for i = 1:nUnitsRest
        for d = 1:max(density)
            if i <= nUnitsRest/2
                biVal(d) = binornd(1,1-biProb);
            else
                biVal(d) = binornd(1,biProb);
            end
        end
        biValInt = biVal*(nUnitsRest/2);
        indW = randperm(nUnitsRest/2,max(density));
        indxLRest(i,:,lay) = indW+biValInt;
        clear biVal biValInt indW
    end
end

