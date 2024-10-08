function [leftEYE , rightEYE , TTL_sInfo] = getEYErawEpoch(rawTIME, ttlTABLE, trialsOfInt, prevtrial, nexttrial , respmat)

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

% trouble shoot whether table has a cell instead of numeric
if iscell(respTABLE.respTime)

    emptyLoc = cellfun(@(x) isempty(x), respTABLE.respTime , 'UniformOutput',true);
    respTABLE.respTime{emptyLoc} = NaN;
    respTABLE.respTime = cell2mat(respTABLE.respTime);

end

allRESPONSE = respTABLE.respTime - respTABLE.QuesStart;

%% Preallocate
TTL_sInfo = cell(numel(unique(ttlTABLE.trialID(ttlTABLE.trialID ~=0 & ~isnan(ttlTABLE.trialID)))),1);

%% Loop through trials
for tttrialir = 1:length(trialsOfInt)

    %% Grab TTLs for the event
    tmpOtr = trialsOfInt(tttrialir);
    tmpOtrTAB = ttlTABLE(ismember(ttlTABLE.trialID, tmpOtr),:);

    % Check if response is MISSING - is NAN in RESPMAT
    if isnan(tmpOtrTAB.TTLid(5))

        % if tmpOtrTAB == 1;

        % If yes populate entire TRIAL ROW with NANs
        columnNames = {'oT_posit_raw', 'oT_posit_cen', 'oT_posit_sd',...
            'oT_posit_cd','oT_posit_dist','oT_pupilS_raw','oT_pupilS_mean',...
            'oT_pupilS_sd','oT_pupilS_cd'};  

        % Create an array of NaNs of the appropriate size
        cellData1 = cell(1, 1);
        cellData2 = cell(1, 1);
        cellData3 = cell(1, 1);
        numericData1 = nan(1, 1);
        cellData4 = cell(1, 1);
        cellData5 = cell(1, 1);
        numericData2 = nan(1, 1);
        numericData3 = nan(1, 1);
        numericData4 = nan(1, 1);

        % Create the table with NaN values and specified column names
        nanTAB = table(cellData1, cellData2, cellData3, numericData1,...
            cellData4, cellData5, numericData2, numericData3,...
            numericData4, 'VariableNames', columnNames);

        % Add to table
        leftEYE(tttrialir,:) = nanTAB;
        rightEYE(tttrialir,:) = nanTAB;

    else


        %% Get the trial times
        trialTimes = double(tmpOtrTAB{:,3});
        trialTTLids = tmpOtrTAB.TTLid;
        if trialTTLids(1) > 100

            % remove the current response row
            tmpOtrTAB = tmpOtrTAB(2:end,:);
            % repopulate response row

            % get info from Respmat
            tmpRESPONSE = respTABLE(tmpOtrTAB.trialID(1),:);
            tmpOtrTAB.ELmessage{5} = ['TTL=',num2str(tmpRESPONSE.respValue{1})];
            tmpOtrTAB.trialID(5) = tmpOtrTAB.trialID(1);
            tmpOtrTAB.timeStamp(5) = round(allRESPONSE(tmpOtrTAB.trialID(1))*1000) + tmpOtrTAB.timeStamp(4);
            tmpOtrTAB.TTLid(5) = tmpRESPONSE.respValue{1};

        end

        diffOFFsets = diff(tmpOtrTAB.timeStamp);
        if diffOFFsets(end) < 20

            tmpOtrTAB.timeStamp(5) = round(allRESPONSE(tmpOtrTAB.trialID(1))*1000) + tmpOtrTAB.timeStamp(4);

        end

        %% Grab next and previous trial event times
        if tmpOtr == 1
            tmpOtrTAB2 = ttlTABLE(ismember(ttlTABLE.trialID, tmpOtr+1),:);
            nexttrialTimes = double(tmpOtrTAB2{:,3});
        elseif tmpOtr > 1 && tmpOtr < length(trialsOfInt)
            tmpOtrTAB2 = ttlTABLE(ismember(ttlTABLE.trialID, tmpOtr+1),:);
            tmpOtrTAB3 = ttlTABLE(ismember(ttlTABLE.trialID, tmpOtr-1),:);
            nexttrialTimes = double(tmpOtrTAB2{:,3});
            prevtrialTimes = double(tmpOtrTAB3{:,3});
        elseif tmpOtr == length(trialsOfInt)
            tmpOtrTAB3 = ttlTABLE(ismember(ttlTABLE.trialID, tmpOtr-1),:);
            prevtrialTimes = double(tmpOtrTAB3{:,3});
        end

        %% Print warnings for trial times and grab what I can
        if exist('nexttrialTimes','var') == 1 && exist('prevtrialTimes','var') == 0
            if nexttrialTimes(1,1) - trialTimes(5,1) <= nexttrial
                numSampnext = nexttrialTimes(1,1) - trialTimes(5,1) - 1;
                warning('Eye tracking interval for trial %s overlaps with next trial. The maximum number of samples you can use is %s.',num2str(tmpOtr),num2str(numSampnext));
                fprintf('\n');
            else
                numSampnext = trialTimes(5,1) + nexttrial;
            end
            startTS = trialTimes(1,1) - 20;
            endTS = trialTimes(5,1) + numSampnext;
            clear nexttrialTimes trialTimes
        elseif exist('nexttrialTimes','var') == 1 && exist('prevtrialTimes','var') == 1
            if prevtrialTimes(5,1) ~= 0
                if trialTimes(1,1) - prevtrialTimes(5,1) <= prevtrial
                    numSampprev = trialTimes(1,1) - prevtrialTimes(5,1) - 1;
                    warning('Eye tracking interval for trial %s overlaps with previous trial. The maximum number of samples you can use is %s.',num2str(tmpOtr),num2str(numSampprev));
                    fprintf('\n');
                else
                    numSampprev = trialTimes(1,1) - prevtrial;
                end
            else
                fprintf('Trial %d is a bad trial..\n',tttrialir-1);
                numSampprev = 5; % Take 5 samples prior to trial onset %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %continue
            end

            if trialTimes(5,1) ~= 0
                if nexttrialTimes(1,1) - trialTimes(5,1) <= nexttrial
                    numSampnext = nexttrialTimes(1,1) - trialTimes(5,1) - 1;
                    warning('Eye tracking interval for trial %s overlaps with next trial. The maximum number of samples you can use is %s.',num2str(tmpOtr),num2str(numSampnext));
                    fprintf('\n');
                else
                    numSampnext = trialTimes(5,1) + nexttrial;
                end
            else
                fprintf('Trial %d is a bad trial..\n',tttrialir);
                numSampnext = 5; % Take 5 samples after trial offset %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %continue
            end
            startTS = trialTimes(1,1) - numSampprev;
            endTS = trialTimes(5,1) + numSampnext;
            clear nexttrialTimes trialTimes prevtrialTimes
        elseif exist('nexttrialTimes','var') == 0 && exist('prevtrialTimes','var') == 1
            if trialTimes(1,1) - prevtrialTimes(5,1) <= prevtrial
                numSampprev = trialTimes(1,1) - prevtrialTimes(5,1) - 1;
                warning('Eye tracking interval for trial %s overlaps with previous trial. The maximum number of samples you can use is %s.',num2str(tmpOtr),num2str(numSampprev));
                fprintf('\n');
            else
                numSampprev = trialTimes(1,1) - prevtrial;
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

        % % Check location of Response - % JAT 9/14/2024
        % if any(diff(tmpTable.ELNKttl) == 0)
        %     % 1. Check order of EventDesc
        %     % 2. Update based on respMat
        %     %
        %     trialROW = respTABLE(tttrialir,:);
        %     durationResp = trialROW.respTime - trialROW.QuesStart; % in seconds
        %     resp2fixcross = trialROW.respTime - trialROW.CrossStart; % in seconds
        %     respDuraSamps = round(durationResp * 1000);
        %     respFixCrSamps = round(resp2fixcross * 1000);
        %
        %     keyboard
        % end


        %% Get event sample number

        disp(tttrialir)
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

    disp(['Finished TRIAL # ',num2str(tttrialir)])

end
end