function [leftEYE , rightEYE , TTL_sInfo] = getEYErawEpoch(rawTIME, ttlTABLE, trialsOfInt, prevtrial, nexttrial)

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

%% Preallocate
TTL_sInfo = cell(numel(unique(ttlTABLE.TrialID(ttlTABLE.TrialID ~=0))),1);

%% Loop through trials
for tttrialir = 1:length(trialsOfInt)

    %% Grab TTLs for the event
    tmpOtr = trialsOfInt(tttrialir);
    tmpOtrTAB = ttlTABLE(ismember(ttlTABLE.TrialID, tmpOtr),:);

    %% Get the trial times
    trialTimes = double(tmpOtrTAB{:,3});

    %% Grab next and previous trial event times
    if tmpOtr == 1
        tmpOtrTAB2 = ttlTABLE(ismember(ttlTABLE.TrialID, tmpOtr+1),:);
        nexttrialTimes = double(tmpOtrTAB2{:,3});
    elseif tmpOtr > 1 && tmpOtr < length(trialsOfInt)
        tmpOtrTAB2 = ttlTABLE(ismember(ttlTABLE.TrialID, tmpOtr+1),:);
        tmpOtrTAB3 = ttlTABLE(ismember(ttlTABLE.TrialID, tmpOtr-1),:);
        nexttrialTimes = double(tmpOtrTAB2{:,3});
        prevtrialTimes = double(tmpOtrTAB3{:,3});
    elseif tmpOtr == length(trialsOfInt)
        tmpOtrTAB3 = ttlTABLE(ismember(ttlTABLE.TrialID, tmpOtr-1),:);
        prevtrialTimes = double(tmpOtrTAB3{:,3});
    end

    %% Print warnings for trial times and grab what I can
    if exist('nexttrialTimes','var') == 1 && exist('prevtrialTimes','var') == 0
        if nexttrialTimes(1,1) - trialTimes(5,1) <= nexttrial
            numSampnext = nexttrialTimes(1,1) - trialTimes(5,1) - 1;
            warning('Eye tracking interval for trial %s overlaps with next trial. The maximum number of samples you can use is %s.',num2str(tmpOtr),num2str(numSampnext));
            fprintf('\n');
        end
        startTS = trialTimes(1,1) - 20;
        endTS = trialTimes(5,1) + numSampnext;
        clear nexttrialTimes trialTimes
    elseif exist('nexttrialTimes','var') == 1 && exist('prevtrialTimes','var') == 1
        if trialTimes(1,1) - prevtrialTimes(5,1) <= prevtrial
            numSampprev = trialTimes(1,1) - prevtrialTimes(5,1) - 1;
            warning('Eye tracking interval for trial %s overlaps with previous trial. The maximum number of samples you can use is %s.',num2str(tmpOtr),num2str(numSampprev));
            fprintf('\n');
        end
        if nexttrialTimes(1,1) - trialTimes(5,1) <= nexttrial
            numSampnext = nexttrialTimes(1,1) - trialTimes(5,1) - 1;
            warning('Eye tracking interval for trial %s overlaps with next trial. The maximum number of samples you can use is %s.',num2str(tmpOtr),num2str(numSampnext));
            fprintf('\n');
        end
        startTS = trialTimes(1,1) - numSampprev;
        endTS = trialTimes(5,1) + numSampnext;
        clear nexttrialTimes trialTimes prevtrialTimes
    elseif exist('nexttrialTimes','var') == 0 && exist('prevtrialTimes','var') == 1
        if trialTimes(1,1) - prevtrialTimes(5,1) <= prevtrial
            numSampprev = trialTimes(1,1) - prevtrialTimes(5,1) - 1;
            warning('Eye tracking interval for trial %s overlaps with previous trial. The maximum number of samples you can use is %s.',num2str(tmpOtr),num2str(numSampprev));
            fprintf('\n');
        end
        startTS = trialTimes(1,1) - numSampprev;
        endTS = trialTimes(5,1) + 20;
        clear prevtrialTimes trialTimes
    end

    %% Create table
    tmpTable = table;
    tmpTable.TrialNUM = repmat(tmpOtr,height(tmpOtrTAB) + 2,1);
    tmpTable.EventID = [nan ; tmpOtrTAB.TTLid ; nan];
    if height(tmpTable.EventID) == 7
        tmpTable.EventDesc = {'Pre';'FixCross'; 'ClipOn';'ClipOff'; 'Probe';'Response';'Post'};
    end
    tmpTable.NLXttl = zeros(height(tmpOtrTAB.TTLid) + 2,1,'int64');
    tmpTable.ELNKttl = [startTS ; tmpOtrTAB.timeStamp ; endTS];
    [tsBlk_OUT] = getTSBlock(startTS,endTS,rawTIME);

    %% Get event sample number
    middleEvents = zeros(height(tmpOtrTAB),1,'int64');
    for mi = 1:height(tmpOtrTAB)
        [~,elnkLoc] = min(abs(double(tmpOtrTAB.timeStamp(mi)) - tsBlk_OUT.Time));
        middleEvents(mi) = elnkLoc;
    end

    %% Add sample number
    tmpTable.ELNKint = [1 ; middleEvents ; height(tsBlk_OUT)];

    %% Move to cell
    TTL_sInfo{tttrialir} = tmpTable;

    %% Grab pupil info
    pupilS_1 = tsBlk_OUT.PupilS_0(:,1);
    pupilS_2 = tsBlk_OUT.PupilS_1(:,1);

    %% Clean up pupil info
    pos_1 = [tsBlk_OUT.GX_0(:,1) , tsBlk_OUT.GY_0(:,1)];
    pos_1c = cleanUPpos(pos_1);
    pos_2 = [tsBlk_OUT.GX_1(:,1) , tsBlk_OUT.GY_1(:,1)];
    pos_2c = cleanUPpos(pos_2);

    %% Create eye table
    [left1 , right2] = createEYEtable(double(pupilS_1),double(pupilS_2),pos_1c,pos_2c);

    %% Add to table
    leftEYE(tttrialir,:) = left1;
    rightEYE(tttrialir,:) = right2;

end
end