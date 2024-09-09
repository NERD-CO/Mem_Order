function respMat = convertrespVals(respMat)

%% function respMat = convertrespVals(respMat)
% This function converts the old confidence values to the new confidence
% values.

%% Set TTL values
origVals = [-3 -2 -1 1 2 3];
newVals = [111 112 113 114 115 116];

%% Loop through respMat
for ii = 1:length(respMat)
    [~,idx] = ismember(respMat(ii).respValue,origVals);
    respMat(ii).respValue = newVals(idx);
end

end