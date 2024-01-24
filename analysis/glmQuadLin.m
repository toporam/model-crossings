function betasAll = glmQuadLin(vects)

if size(vects,2) == 5
    v = [-1 -.5 0 .5 1];
elseif size(vects,2) == 3
    v = [-1 0 1];
end 

r1 = v/norm(v);
y = v.^2; 
y = y - mean(y);
r2 = y/norm(y);

X = [r1', r2'];

betasAll = glmfit(X,vects');