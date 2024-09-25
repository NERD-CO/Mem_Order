function [corValsNew] = computeASCI_eye_v2(waveforms)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

% if numel(unique(idx)) == 1;
%      corInds = nan(length(idx),1);
%      return
% end


% tempWaves = transpose(waveforms);
tempWaves = waveforms;
% corValsOld = zeros(size(tempWaves,1),1);
corValsNew = zeros(size(tempWaves,1),1);

% Create Rp, Rn and Rz from clust template
meanCw = mean(tempWaves,'omitnan');
template = meanCw/max(abs(meanCw));

% get start number of waveforms
startNumWaves = size(tempWaves,1);

% crit = 0.625;

% Start While loop

Rp = nan(size(template));
Rn = nan(size(template));

% for mi = 1:length(template)
% 
%     tRz = template(mi);
% 
%     if abs(tRz) <= 0.20
%         Rp(mi) = 0.20;
%         Rn(mi) = -0.20;
%     elseif tRz > 0.20
%         Rp(mi) = max([tRz/2, 0.20]);
%         Rn(mi) = Rp(mi) - 0.20;
%     else
%         Rn(mi) = min([tRz/2 , -0.20]);
%         Rp(mi) = Rn(mi) + 0.20;
%     end
% 
% end
% 
% Rp2 = Rp + template;
% Rn2 = template - Rn;

M = movmad(template,20);
M2 = M*75;
% plot(template + (M*10)); hold on; plot(template - (M*10))
Rp2 = M2 + template;
Rn2 = template - M2;

for ai = 1:startNumWaves

    tWave = tempWaves(ai,:);
    tWaveN = tWave/max(abs(meanCw));

    SY = nan(1,size(tWaveN,2));
    % SX = template;

    for ti = 1:size(tWaveN,2)

        yi = tWaveN(ti);
        if yi <= Rp2(ti) && yi >= Rn2(ti)
            SY(ti) = 1;
        elseif yi >= Rp2(ti)
            SY(ti) = 0;
        elseif yi <= Rn2(ti)
            SY(ti) = 0;
        elseif isnan(yi)
            SY(ti) = 0;
        end

        % xi = template(ti);
        % if xi <= Rp2(ti) && xi >= Rn2(ti)
        %     SX(ti) = 0;
        % elseif xi >= Rp2(ti)
        %     SX(ti) = 1;
        % elseif xi <= Rn2(ti)
        %     SX(ti) = -1;
        % elseif isnan(xi)
        %     SX(ti) = 0;
        % end

    end

    %------------------------------------------------------------
    % To satisfy dot prod makes indicies where both zeros equal 1
    % bzInd = SY == 0 & SX == 0;
    % SX(bzInd) = 1;
    % SY(bzInd) = 1;
    %------------------------------------------------------------

    %------------------------------------------------------------
    % Factor in NaNs to length
    SYlen = sum(~isnan(tWaveN));

    %------------------------------------------------------------

    % SCC1 = sum(SY.*SX)/length(SY); 
    SCC2 = sum(SY)/SYlen; % was length(SY)

    % corValsOld(ai) = SCC1;
    corValsNew(ai) = SCC2;

end



test = 1;



end

