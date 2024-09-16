function respMat = convertrespVals(respMat)

%% function respMat = convertrespVals(respMat)
% This function converts the old confidence values to the new confidence
% values.

%% Set TTL values
origVals = [-3 -2 -1 1 2 3];
newVals = [111 112 113 114 115 116];

%% Loop through respMat
for ii = 1:length(respMat)
    if ismember(respMat(ii).respValue,origVals)
        [~,idx] = ismember(respMat(ii).respValue,origVals);
        respMat(ii).respValue = newVals(idx);
    else
        respMat(ii).respValue = NaN;
    end
end

end