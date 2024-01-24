function allLayersGF = applyGF(allLayers, allGF)

nLayers = length(allLayers);

allLayersGF{1} = allGF{1} .* allLayers{1}; %multiply GF by the activation values for layer 1 

for i = 2:nLayers 
    allLayersGF{i} = allGF{i} .* allLayers{i};  % same for the remaining layers
end

