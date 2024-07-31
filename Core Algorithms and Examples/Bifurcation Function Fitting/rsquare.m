function out = rsquare(y,ypred)

% Compute r-square of predictions

% Make both row vector
if size(y,1)>1
    y = y';
end
if size(ypred,1)>1
    ypred = ypred';
end

vartotal = var(y)*(length(y)-1);    % Variance sum

varexp = sum((y-ypred).^2);

out = 1-varexp/vartotal;










