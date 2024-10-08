function [leftEYE , rightEYE , TTL_sInfo] = getEYErawEpoch_V2(rawTIME, ttlTABLE, trialsOfInt, prevtrial, nexttrial , respmat)

%% function [leftEYE , rightEYE , TTL_sInfo] = getEYErawEpoch(rawTIME , ttlTABLE , trialsOfInt)
% NEW RAW will go from -400ms [from start][total 500 ITI] to +200ms [from end]
% NEW TABLE per trial: EVENT , TRIALnum , NLX_T [future] , EYELink_T
% TABLE: TRIAL_Num , EVENT_ID , NLX_TTL , EYElink_TTL , EYElink_Int

%% Make N by 2 matrix of fieldname + value type
bothVarNames = [['oT_posit_raw', "cell"]; ...
    ['oT_posit_cen', "cell"]; ...
    ['oT_posit_sd', "cell"]; ...
    ['oT_posit_cd', "double"]; ...
    ['oT_posit_dist', "cell"]; ...
    ['oT_pupilS_raw', "cell"];...
    ['oT_pupilS_mean', "double"];...
    ['oT_pupilS_sd', "double"];...
    ['oT_pupilS_cd', "double"]];

%% Make table using fieldnames & value types from above
leftEYE = table('Size',[0,size(bothVarNames,1)],...
    'VariableNames', bothVarNames(:,1),...
    'VariableTypes', bothVarNames(:,2));

rightEYE = table('Size',[0,size(bothVarNames,1)],...
    'VariableNames', bothVarNames(:,1),...
    'VariableTypes', bothVarNames(:,2));

%% Convert respMat to table
respmat = convertrespVals(respmat);
respTABLE = struct2table(respmat);
allRESPONSE =respTABLE.respTime - respTABLE.QuesStart;

%% Preallocate
TTL_sInfo = cell(length(respmat),1);

%% Loop to calculate time differences in respMat
for ii = 1:length(respmat)
    crossTime(ii,1) = respmat(ii).CrossEnd - respmat(ii).CrossStart;
    crossITI(ii,1) = respmat(ii).FrameStart - respmat(ii).CrossEnd;
    frameTime(ii,1) = respmat(ii).FrameEnd - respmat(ii).FrameStart;
    frameITI(ii,1) = respmat(ii).QuesStart - respmat(ii).FrameEnd;
    respTime(ii,1) = respmat(ii).respTime - respmat(ii).QuesStart;
    respValue(ii,1) = respmat(ii).respValue;
end

%% Calculate difference in ttlTable
timeStamp = ttlTABLE{:,3};
timeDiff = diff(timeStamp);
end