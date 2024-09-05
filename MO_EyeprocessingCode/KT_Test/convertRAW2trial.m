function [outTABLE] = convertRAW2trial(inTTLtab)

%% function [outTABLE] = convertRAW2trial(inTTLtab)

%% TTL info
ttlINFO = inTTLtab.taskinformation;
indexLOG = matches(ttlINFO.TTLvalue,'11'); % changed from 1 to 11 - KT 08/27/24
indexLOC = find(indexLOG);
trialCell = cell(sum(indexLOG),1);
trialInVec = zeros(sum(~matches(ttlINFO.TTLvalue,'66')),1);

%% Loop over indices
for ii = 1:numel(indexLOC)
    if ii ~= length(indexLOC)
        tmpStart = indexLOC(ii);
        tmpEnd = indexLOC(ii + 1) - 1;
        trialInVec(tmpStart:tmpEnd) = ii;
        tmpTAB = ttlINFO(tmpStart:tmpEnd,:);
        tmpTAB.TrialID = trialInVec(tmpStart:tmpEnd);
        trialCell{ii} = tmpTAB;
    else
        tmpStart = indexLOC(ii);
        tmpEnd = sum(~matches(ttlINFO.TTLvalue,'66'));
        trialInVec(tmpStart:tmpEnd) = ii;
        tmpTAB = ttlINFO(tmpStart:tmpEnd,:);
        tmpTAB.TrialID = trialInVec(tmpStart:tmpEnd);
        trialCell{ii} = tmpTAB;
    end
end

%% Create out table
outTABLE = trialCell;

end