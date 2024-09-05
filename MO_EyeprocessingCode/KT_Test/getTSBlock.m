function [tsBlk_OUT] = getTSBlock(startI,endI,rawT)

%% function [tsBlk_OUT] = getTSBlock(startI,endI,rawT)

%% Get ts info
[~, eyeTTL1_i] = min(abs(double(startI) - rawT.Time));
[~, eyeTTL2_i] = min(abs(double(endI) - rawT.Time));

tsBlk_OUT = rawT(eyeTTL1_i:eyeTTL2_i,:);

end