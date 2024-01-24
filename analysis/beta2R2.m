function [varExp] = beta2R2(betas, vect, X)

if length(betas)>3
    warning('input vector largestDim > 3, which should not happen'),
end

varExp(1) = var(X(:,1).*betas(2))/var(vect);
varExp(2) = var(X(:,2).*betas(3))/var(vect);



