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

%% Fix tsTable
tsTable = Fix_MemOrder_tsTable(tsTable,behavioR);

%% Modify tsTable to remove extra TTLs
responseEvents = [111 112 113 114 115 116];
% tossEntry = nan(length(tsTable.TTLid),1);
% for ii = 1:length(tsTable.TTLid)
%     if ~ismember(responseEvents,tsTable.TTLid(ii,1))
%         tossEntry(ii,1) = 0;
%     else
%         num1 = tsTable.TTLid(ii-1,1);
%         num2 = tsTable.TTLid(ii,1);
%         if ismember(num1,responseEvents) && ismember(num2,responseEvents)
%             tossEntry(ii,1) = 1;
%         else
%             tossEntry(ii,1) = 0;
%         end
%     end
% end
% tossEntry = logical(tossEntry);
% tsTable = tsTable(~tossEntry,:);

%% Add trial ID
% taskEVEnts = [11 1 2 3 4 0];
% trialIDtmp = zeros(height(tsTable),1);
% for ttIi = 1:numel(taskEVEnts)
%     TTLIDlocs = find(tsTable.TTLid == taskEVEnts(ttIi));
%     TTLnuMS = 1:numel(TTLIDlocs);
%     trialIDtmp(TTLIDlocs) = TTLnuMS;
% end
% 
% responseEvents = [111 112 113 114 115 116];
% [eventLoc,~] = ismember(tsTable.TTLid,responseEvents);
% TTLnuMS = 1:sum(eventLoc);
% trialIDtmp(eventLoc) = TTLnuMS;
% 
% tsTable.TrialID = trialIDtmp;

%% Add rows for responses if they don't exist
if ~any(ismember(tsTable.TTLid,responseEvents))
    response = [-3 -2 -1 1 2 3];
    newRow = table(NaN,NaN,NaN,NaN,'VariableNames',tsTable.Properties.VariableNames);
    idx = find(tsTable.TTLid == 3);
    for i = length(idx):-1:1

        %% Add empty row
        tsTable = [tsTable(1:idx(i),:); newRow; tsTable(idx(i)+1:end,:)];

        %% Create ELmessage
        [~,val] = ismember(behavioR(i).respValue,response);
        tsTable.ELmessage{idx(i)+1} = append('TTL=',num2str(responseEvents(val)));

        %% Add TTL
        if ~isempty(val)
            tsTable.TTLid(idx(i)+1) = responseEvents(val);
        end

        %% Add trial ID
        tsTable.TrialID(idx(i)+1) = tsTable.TrialID(idx(i));

        %% Add time stamp
        eventTime = (behavioR(i).respTime - behavioR(i).QuesStart)*1000;
        tsTable.timeStamp(idx(i)+1) = tsTable.timeStamp(idx(i)) + eventTime;

        %% Warn if time issue
        if tsTable.timeStamp(idx(i)+1) > tsTable.timeStamp(idx(i)+2)
            warning('Eye tracking for trial %s overlaps with trial %s..\n',i,i+1)
        end

    end
end

%% Create behavior table
clipStime = tsTable.timeStamp(tsTable.TTLid == 2);
clipStrial = tsTable.TrialID(tsTable.TTLid == 2);
if isfield(behavioR,'ClipName') % Encode
    clipIDS = {behavioR.ClipName};
elseif isfield(behavioR,'FrameName') % Scene Recog
    clipIDS = {behavioR.FrameName};
end
quesStime = tsTable.timeStamp(tsTable.TTLid == 3);
if isfield(behavioR,'QuesName') % Encode
    quesIDS = {behavioR.QuesName};
else
    quesIDS = cell(1,length(behavioR));
end
behTABLE = struct2table(behavioR);
responTimeSecs = behTABLE.respTime - behTABLE.QuesStart;

%% Create vidQuestable
vidQuestable = table(clipStrial ,clipStime,transpose(clipIDS), quesStime ,...
    transpose(quesIDS),...
    responTimeSecs,'VariableNames',{'TrialNum','ClipTstamp','ClipName',...
    'QuestTstamp','QuestName','ResponseSecs'});

end