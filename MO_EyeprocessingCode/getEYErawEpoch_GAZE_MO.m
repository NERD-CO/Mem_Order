function [leftEYE , rightEYE , TTL_sInfo] =...
    getEYErawEpoch_GAZE_MO(rawTIME, ttlTABLE, trialsOfInt,...
    picTABLE , picLOCation, prevtrial, nexttrial , respmat)


%% Make N by 2 matrix of fieldname + value type
bothVarNames = [['gaze_Raw', "cell"]; ...
    ['gaze_Raw_cen', "cell"]; ...
    ['gaze_Raw_sd', "cell"]; ...
    ['gaze_Scr', "cell"]; ...
    ['gaze_Scr_cen', "cell"];...
    ['gaze_Scr_sd', "cell"];...
    ['gaze_Pic', "cell"];...
    ['gaze_Pic_cen', "cell"];...
    ['gaze_Pic_sd', "cell"];...
    ['gaze_Scr_time', "cell"];
    ['gaze_Pic_time',"cell"]];

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
        columnNames = {'gaze_Raw', 'gaze_Raw_cen', 'gaze_Raw_sd',...
            'gaze_Scr','gaze_Scr_cen','gaze_Scr_sd','gaze_Pic',...
            'gaze_Pic_cen','gaze_Pic_sd','gaze_Scr_time','gaze_Pic_time'};  

        % Create an array of NaNs of the appropriate size
        cellData1 = cell(1, 1);
        cellData2 = cell(1, 1);
        cellData3 = cell(1, 1);
        numericData1 = cell(1, 1);
        cellData4 = cell(1, 1);
        cellData5 = cell(1, 1);
        numericData2 = cell(1, 1);
        numericData3 = cell(1, 1);
        numericData4 = cell(1, 1);
        cellData6 = cell(1, 1);
        cellData7 = cell(1, 1);

        % Create the table with NaN values and specified column names
        nanTAB = table(cellData1, cellData2, cellData3, numericData1,...
            cellData4, cellData5, numericData2, numericData3,...
            numericData4,cellData6,cellData7, 'VariableNames', columnNames);

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
                continue
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
                continue
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

        %% load in first image of screen / image / or pairs of images
        % May need to separate out ENcoding, SCeneRe, 
        picID = [picLOCation , filesep , picTABLE.ClipName{tttrialir}];

        [xPictBnds , yPictBnds] = getPICtureDims_MO(picLOCation ,  picTABLE.ClipName{tttrialir});


        %% Create eye table
        [left1 , right2] = createEYEtable(double(pupilS_1),double(pupilS_2),pos_1c,pos_2c);

        %% Add to table
        leftEYE(tttrialir,:) = left1;
        rightEYE(tttrialir,:) = right2;

    end

    disp(['Finished TRIAL # ',num2str(tttrialir)])

end
end