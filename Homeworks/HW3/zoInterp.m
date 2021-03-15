function [copies] = zoInterp(x, numInterp)

copies = reshape(repmat(x, numInterp, 1), 1, length(x) * numInterp);

end

