function [x1,x2] = patchNaN(x1,x2)

% This function is used to patch NaNs to make x1 and x2 of the same length

% Column vector is needed

if length(x1)>length(x2)
    x2 = [x2; NaN(length(x1)-length(x2),1)];
elseif length(x1)<length(x2)
    x1 = [x1; NaN(length(x2)-length(x1),1)];
end