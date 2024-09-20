function [tsTable2] = Fix_MemOrder_tsTable(tsTable,behavioR)

%% function [tsTable] = Fix_MemOrder_tsTable(tsTable)
% Notes for later

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

    if ismember(tmpTTLid,[60,61])
        trialID(i) = nan;
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

uniTRAILid = unique(cleanTABLEt.trialID(cleanTABLEt.trialID ~= 0 & ~isnan(cleanTABLEt.trialID)));
behavioR2 = behavioR;

behavioR2 = convertrespVals(behavioR2);
respMATT = struct2table(behavioR2);

% Check respTime if cell or not - added 9/20/2024
if iscell(respMATT.respTime)

    emptyLoc = cellfun(@(x) isempty(x), respMATT.respTime , 'UniformOutput',true);
    respMATT.respTime{emptyLoc} = NaN;
    respMATT.respTime = cell2mat(respMATT.respTime);

end

chkRESPONSE = zeros(numel(uniTRAILid),4);
for uui = 1:height(respMATT)

    tmpRESPONSE = respMATT.respValue(uui);
    tmpResponseT = respMATT.respTime(uui) - respMATT.QuesStart(uui);
    tmpResponseS = round(tmpResponseT*1000);

    if isnan(tmpRESPONSE)
        continue
    else

        c_trialTable = cleanTABLEt(ismember(cleanTABLEt.trialID,uui),:);
        if uui == length(uniTRAILid)

            startINDEX = c_trialTable.startTimesA(4) + 1;
            endINDEX = height(tsTable) - 1;
        else
            n_trialTable = cleanTABLEt(ismember(cleanTABLEt.trialID,uui + 1),:);

            startINDEX = c_trialTable.startTimesA(4) + 1;
            endINDEX = n_trialTable.startTimesA(1) - 1;
        end


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

    end

end

%% Add responses for tsTable
tsTable2 = cleanTABLEt;
qIdx = find(tsTable2.TTLid == 3);
newRow = table(NaN,NaN,NaN,NaN,NaN,'VariableNames',cleanTABLEt.Properties.VariableNames);
for i = length(qIdx):-1:1

    %% Add empty row
    tsTable2 = [tsTable2(1:qIdx(i),:); newRow; tsTable2(qIdx(i)+1:end,:)];

    %% Create ELmessage
    val = chkRESPONSE(i,1);
    if val == 0
        tsTable2.ELmessage{qIdx(i)+1} = append('TTL=','NaN');
    else
        tsTable2.ELmessage{qIdx(i)+1} = append('TTL=',num2str(val));
    end

    %% Add TTL
    if val == 0
        tsTable2.TTLid(qIdx(i)+1) = NaN;
    else
        tsTable2.TTLid(qIdx(i)+1) = val;
    end

    %% Add trial ID
    tsTable2.trialID(qIdx(i)+1) = tsTable2.trialID(qIdx(i));

    %% Add time stamp
    if val == 0
        tsTable2.timeStamp(qIdx(i)+1) = NaN;
    else
        tsTable2.timeStamp(qIdx(i)+1) = chkRESPONSE(i,3);
    end

    %% Add NaN for start time
    tsTable2.startTimesA(qIdx(i)+1) = NaN;
    
end
tsTable2.startTimesA = [];

end