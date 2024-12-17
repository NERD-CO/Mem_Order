% tsT_top = tsTable.timeStamp(1);
% tsT_bot = tsTable.timeStamp(end);
% 
% duratiONN = tsT_bot - tsT_top;
% 
% totalTime = duratiONN/1000/60;

%%
% 1. Extract indices for 11
% 2. Extract time points for 11
% 3. Compute offset between timestamps 
% 4. Compute offset between TTLs
% 5. Compute number of samples between intervening TTLs
% 6. Prep for individual trial

% % 1
% startIndicies = find(ismember(tsTable.TTLid,11));
% % 2
% startTimes = tsTable.timeStamp(startIndicies);
% % 3
% diff(double(startTimes))/1000;
% % 4
% diff(startIndicies);

%%
startIndicies = find(ismember(tsTable.TTLid,11));
questIndicies = find(ismember(tsTable.TTLid,3));
cleanTABLEt = tsTable(ismember(tsTable.TTLid,[11 1 2 3 61 60]),:);

startTimesA = nan(height(cleanTABLEt),1);

startTimesA(cleanTABLEt.TTLid == 11) = startIndicies;
startTimesA(cleanTABLEt.TTLid == 3) = questIndicies;

cleanTABLEt.startTimesA = startTimesA;

% Create trial IDS
trialID = nan(height(cleanTABLEt),1);
curTrial = 0;
for i = 1:height(cleanTABLEt)

    tmpTTLid = cleanTABLEt.TTLid(i);

    if ismember([60 61],tmpTTLid)
        continue
    end

    if tmpTTLid == 11
        curTrial = curTrial + 1;
        trialID(i) = curTrial;
    else
        trialID(i) = curTrial;
    end

end

cleanTABLEt.trialID = trialID;

%% Explore response data

uniTRAILid = unique(cleanTABLEt.trialID(cleanTABLEt.trialID ~= 0));
behavioR2 = behavioR;

behavioR2 = convertrespVals(behavioR2);
respMATT = struct2table(behavioR2);

perRESPONSE = zeros(numel(uniTRAILid),1);
chkRESPONSE = zeros(numel(uniTRAILid),4);
for uui = 1:height(respMATT)

    tmpRESPONSE = respMATT.respValue(uui);
    tmpResponseT = respMATT.respTime(uui) - respMATT.QuesStart(uui);
    tmpResponseS = round(tmpResponseT*1000);

    if isnan(tmpRESPONSE)
        perRESPONSE(uui) = nan;
    else

        c_trialTable = cleanTABLEt(ismember(cleanTABLEt.trialID,uui),:);
        if uui == length(uniTRAILid)
            n_trialTable = cleanTABLEt(ismember(cleanTABLEt.trialID,height(respMATT)),:);
        else
            n_trialTable = cleanTABLEt(ismember(cleanTABLEt.trialID,uui + 1),:);
        end
        startINDEX = c_trialTable.startTimesA(4) + 1;
        endINDEX = n_trialTable.startTimesA(1) - 1;

        % Get intervening responses
        interTable = tsTable(startINDEX:endINDEX,:);

        % Sample duration
        sampleDurF = double(c_trialTable.timeStamp(4) + tmpResponseS);

        % Find indices of correct responses
        % interTableU = interTable(interTable.TTLid == tmpRESPONSE,:);

        % Find indices near expected range
        [~ , mloc] = min(abs(double(interTable.timeStamp) - sampleDurF));
        tmpCHECK = interTable.TTLid(mloc);

        chkRESPONSE(uui,1) = tmpCHECK;
        chkRESPONSE(uui,2) = tmpRESPONSE;
        chkRESPONSE(uui,3) = interTable.timeStamp(mloc);
        chkRESPONSE(uui,4) = interTable.timeStamp(mloc) - c_trialTable.timeStamp(4);


        responsTABu = cell2table(tabulate(categorical(interTable.TTLid)),...
            'VariableNames',{'respVal','count','percent'});

        responsTABu.respVal = cellfun(@(x) str2double(x) , responsTABu.respVal,...
            'UniformOutput',true);

        perRESPONSE(uui) = responsTABu.percent(ismember(responsTABu.respVal,tmpRESPONSE));
    end

end

perRESPONSE = round(perRESPONSE);
histogram(perRESPONSE)

%%
% 11 cross onset 
% 3 question onset 
% 1. Get repmat response time

% 2. Obtain timestamps from matching responses

% 3. ?

