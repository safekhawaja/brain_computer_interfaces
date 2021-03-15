function [features] = MovingWinFeats(x, fs, winLen, winDisp, featF)

NumWindows = floor((((length(x) / fs) - winLen)/winDisp)+1);
strt = 1;
features = zeros(1, NumWindows);

for i = 1:NumWindows
    ib = x(strt:strt + winLen * fs - 1);
    features(i) = featF(ib);
    strt = strt + winDisp * fs;
end

end

