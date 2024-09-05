function [tsTable,vidQuestable] = createTStable(edfR , behavioR)

%%
mseCol = {edfR.FEVENT.message};
tsCol = [edfR.FEVENT.sttime];
mseCol_Char = cellfun(@(x) char(x), mseCol, "UniformOutput",false);

%% ts table
% Timestamp, TTL ID, raw message
mseColTT = contains(mseCol_Char,'TTL');
messTT = mseCol(mseColTT);
timeStamp1 = tsCol(mseColTT);
ttlNum1 = extractAfter(messTT,"=");
ttlNum2 = cellfun(@(x) str2double(x), ttlNum1, "UniformOutput", true);

tsTable = table(transpose(messTT),transpose(ttlNum2),transpose(timeStamp1),...
    'VariableNames',{'ELmessage','TTLid','timeStamp'});

%% Add trial ID
taskEVEnts = [11 1 2 3 4 0];
trialIDtmp = zeros(height(tsTable),1);
for ttIi = 1:numel(taskEVEnts)
    TTLIDlocs = find(tsTable.TTLid == taskEVEnts(ttIi));
    TTLnuMS = 1:numel(TTLIDlocs);
    trialIDtmp(TTLIDlocs) = TTLnuMS;
end

responseEvents = [111 112 113 114 115 116];
[eventLoc,~] = ismember(tsTable.TTLid,responseEvents);
TTLnuMS = 1:sum(eventLoc);
trialIDtmp(eventLoc) = TTLnuMS;

tsTable.TrialID = trialIDtmp;

%% Create behavior table
clipStime = tsTable.timeStamp(tsTable.TTLid == 2);
clipStrial = tsTable.TrialID(tsTable.TTLid == 2);
clipIDS = {behavioR.ClipName};
quesStime = tsTable.timeStamp(tsTable.TTLid == 3);
quesIDS = {behavioR.QuesName};
behTABLE = struct2table(behavioR);
responTimeSecs = behTABLE.respTime - behTABLE.QuesStart;

%% Create vidQuestable
vidQuestable = table(clipStrial ,clipStime,transpose(clipIDS), quesStime ,...
    transpose(quesIDS),...
    responTimeSecs,'VariableNames',{'TrialNum','ClipTstamp','ClipName',...
    'QuestTstamp','QuestName','ResponseSecs'});

end