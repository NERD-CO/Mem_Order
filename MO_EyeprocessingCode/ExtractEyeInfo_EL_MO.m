function [tsTable, vidQuestable, fixTABLE_Eye0 ,...
          fixTABLE_Eye1 ,...
          sacTABLE_Eye0 ,...
          sacTABLE_Eye1 ,...
          allgazePosPupsTAB] = ExtractEyeInfo_EL_MO(eyeMatfile , behavOUTPUT)

% ---Function 'ExtractEyeInfo'---
% 
% Extract eyetracking data from converted .edf to .mat file
% 
% Input: eyeMatfile
%   -Generated output from 'Edf2Mat_UCH.m' function
% Output:
%   -tsTable: timestamps (unit = 1000 Hz)
%   -picTable: where picture stimuli are located on behavioral laptop
%       *need to change dir for local machine
%   -fixTab: L and R eye info, time when subject fixates on a certain point
%   -saccTab: rows = unique saccade events
%   -rawTab: clock data
% 
% John A. Thompson | May 29 2022

load(eyeMatfile,'edfRawStruct');
edfRAW = edfRawStruct;

% Picture and Timestamp Tables
[tsTable, vidQuestable] = createTStable(edfRAW , behavOUTPUT);

edfRAWEvent = struct2table(edfRAW.FEVENT);
gxALL = edfRAW.FSAMPLE.gx;
gyALL = edfRAW.FSAMPLE.gy;
timeALL = edfRAW.FSAMPLE.time;
posX = edfRAW.FSAMPLE.px;
posY = edfRAW.FSAMPLE.py;
velX = edfRAW.FSAMPLE.gxvel;
velY = edfRAW.FSAMPLE.gyvel;
pupilS = edfRAW.FSAMPLE.pa;

% allgazeTAB = array2table([transpose(gxALL) , transpose(gyALL) ,...
%     transpose(timeALL)], 'VariableNames', {'GX_0','GX_1',...
%     'GY_0','GY_1','TimeTT'});

allgazePosPupsTAB = array2table([transpose(timeALL) , transpose(posX) , transpose(posY) ,...
    transpose(velX), transpose(velY) , transpose(gxALL),...
    transpose(gyALL) , transpose(pupilS)],'VariableNames',...
    {'Time','PosX_0','PosX_1','PosY_0','PosY_1','VelX_0','VelX_1','VelY_0',...
    'VelY_1','GX_0','GX_1','GY_0','GY_1','PupilS_0','PupilS_1'});

eyeIDs = [0 1];
allEYEs = edfRAWEvent.eye;
eyesOUT = struct;
for eyeI = 1:2

    tmpEYE = eyeIDs(eyeI);
    eyeTAB = edfRAWEvent(allEYEs == tmpEYE,:);

    % Get FIX start
    endFIX = matches(eyeTAB.codestring,'ENDFIX');
    eFIXtab = eyeTAB(endFIX,:);

    eyesOUT.(['eye_',num2str(tmpEYE)]).Fixes = [eFIXtab.sttime,...
        eFIXtab.entime];
end

eyesOUT.eye_0.Fixes = array2table(eyesOUT.eye_0.Fixes,'VariableNames',...
    {'StartTime','EndTime'});
eyesOUT.eye_1.Fixes = array2table(eyesOUT.eye_1.Fixes,'VariableNames',...
    {'StartTime','EndTime'});


for eyeI = 1:2

    tmpEYE = eyeIDs(eyeI);
    eyeTAB = edfRAWEvent(allEYEs == tmpEYE,:);

    % Get FIX start
    endSAC = matches(eyeTAB.codestring,'ENDSACC');
    eSACtab = eyeTAB(endSAC,:);

    eyesOUT.(['eye_',num2str(tmpEYE)]).Saccades = [eSACtab.sttime,...
        eSACtab.entime];
end

eyesOUT.eye_0.Saccades = array2table(eyesOUT.eye_0.Saccades,'VariableNames',...
    {'StartTime','EndTime'});
eyesOUT.eye_1.Saccades = array2table(eyesOUT.eye_1.Saccades,'VariableNames',...
    {'StartTime','EndTime'});


% Process raw gx and gy positions
[fixTABLE_Eye0 ,...
          fixTABLE_Eye1 ,...
          sacTABLE_Eye0 ,...
          sacTABLE_Eye1] = getFIXSACdata(eyesOUT , gxALL , gyALL , timeALL);



end


function [tsTable,vidQuestable] = createTStable(edfR , behavioR)


mseCol = {edfR.FEVENT.message};
tsCol = [edfR.FEVENT.sttime];
mseCol_Char = cellfun(@(x) char(x), mseCol, "UniformOutput",false);

% ts table
% Timestamp, TTL ID, raw message
mseColTT = contains(mseCol_Char,'TTL');
messTT = mseCol(mseColTT);
timeStamp1 = tsCol(mseColTT);
ttlNum1 = extractAfter(messTT,"=");
ttlNum2 = cellfun(@(x) str2double(x), ttlNum1, "UniformOutput", true);

tsTable = table(transpose(messTT),transpose(ttlNum2),transpose(timeStamp1),...
    'VariableNames',{'ELmessage','TTLid','timeStamp'});

taskEVEnts = [11 1 2 3 4 0];

trialIDtmp = zeros(height(tsTable),1);
for ttIi = 1:numel(taskEVEnts)
    TTLIDlocs = find(tsTable.TTLid == taskEVEnts(ttIi));
    TTLnuMS = 1:numel(TTLIDlocs);
    trialIDtmp(TTLIDlocs) = TTLnuMS;
end
tsTable.TrialID = trialIDtmp;

% trialID = unique(tsTable.TrialID(tsTable.TrialID ~= 0));
% CLIP ONSET = 2
% PROBE ONSET = 3
% % picture table
% % Timestamp, picture name, Trial number, 

clipStime = tsTable.timeStamp(tsTable.TTLid == 2);
clipStrial = tsTable.TrialID(tsTable.TTLid == 2);
clipIDS = {behavioR.ClipName};

quesStime = tsTable.timeStamp(tsTable.TTLid == 3);
% quesStrial = tsTable.TrialID(tsTable.TTLid == 3);
quesIDS = {behavioR.QuesName};

behTABLE = struct2table(behavioR);

responTimeSecs = behTABLE.respTime - behTABLE.QuesStart;

vidQuestable = table(clipStrial ,clipStime,transpose(clipIDS), quesStime ,...
    transpose(quesIDS),...
    responTimeSecs,'VariableNames',{'TrialNum','ClipTstamp','ClipName',...
    'QuestTstamp','QuestName','ResponseSecs'});
end











function [fixTABLE_Eye0 ,...
          fixTABLE_Eye1 ,...
          sacTABLE_Eye0 ,...
          sacTABLE_Eye1] = getFIXSACdata(eyesOUT , gxALL , gyALL , timeALL)

% Fixation extract

eyeIDs = [0 1];
fix_eyeZEROp = cell(height(eyesOUT.eye_0.Fixes),1);
fix_eyeZEROd = nan(height(eyesOUT.eye_0.Fixes),1);
fix_eyeONEp = cell(height(eyesOUT.eye_1.Fixes),1);
fix_eyeONEd = nan(height(eyesOUT.eye_1.Fixes),1);
for eye2u = 1:2
    tmpEYE = eyeIDs(eye2u);

    tmpFIXES = eyesOUT.(['eye_',num2str(tmpEYE)]).Fixes;

    for fixI = 1:height(tmpFIXES)
        startINDf = tmpFIXES.StartTime(fixI);
        endINDf = tmpFIXES.EndTime(fixI);

        startTTf = find(timeALL == startINDf);
        endTTf = find(timeALL == endINDf);

        switch tmpEYE
            case 0
                fix_eyeZEROp{fixI} =...
                    array2table([transpose(gxALL(1,startTTf:endTTf))...
                    transpose(gyALL(1,startTTf:endTTf))...
                    transpose(timeALL(1,startTTf:endTTf))],...
                    'VariableNames',{'gx','gy','timeTTL'});

                fix_eyeZEROd(fixI) = (endTTf - startTTf)/1000;

            case 1
                fix_eyeONEp{fixI} =...
                    array2table([transpose(gxALL(1,startTTf:endTTf))...
                    transpose(gyALL(1,startTTf:endTTf))...
                    transpose(timeALL(1,startTTf:endTTf))],...
                    'VariableNames',{'gx','gy','timeTTL'});

                fix_eyeONEd(fixI) = (endTTf - startTTf)/1000;

        end
    end
end

fixTABLE_Eye0 = table(fix_eyeZEROp , fix_eyeZEROd,...
    'VariableNames',{'FIX_Eye0_Points','FIX_Eye0_DurSecs'});

fixTABLE_Eye1 = table(fix_eyeONEp , fix_eyeONEd,...
    'VariableNames',{'FIX_Eye1_Points','FIX_Eye1_DurSecs'});

sac_eyeZEROp = cell(height(eyesOUT.eye_0.Saccades),1);
sac_eyeZEROd = nan(height(eyesOUT.eye_0.Saccades),1);
sac_eyeONEp = cell(height(eyesOUT.eye_1.Saccades),1);
sac_eyeONEd = nan(height(eyesOUT.eye_1.Saccades),1);
for eye2u = 1:2
    tmpEYE = eyeIDs(eye2u);

    tmpSACCADES = eyesOUT.(['eye_',num2str(tmpEYE)]).Saccades;

    for sacI = 1:height(tmpSACCADES)
        startINDs = tmpSACCADES.StartTime(sacI);
        endINDs = tmpSACCADES.EndTime(sacI);

        startTTs = find(timeALL == startINDs);
        endTTs = find(timeALL == endINDs);

        switch tmpEYE
            case 0
                sac_eyeZEROp{sacI} =...
                    array2table([transpose(gxALL(1,startTTs:endTTs))...
                    transpose(gyALL(1,startTTs:endTTs))...
                    transpose(timeALL(1,startTTs:endTTs))],...
                    'VariableNames',{'gx','gy','timeTTL'});

                sac_eyeZEROd(sacI) = (endTTs - startTTs)/1000;

            case 1
                sac_eyeONEp{sacI} =...
                    array2table([transpose(gxALL(1,startTTs:endTTs))...
                    transpose(gyALL(1,startTTs:endTTs))...
                    transpose(timeALL(1,startTTs:endTTs))],...
                    'VariableNames',{'gx','gy','timeTTL'});

                sac_eyeONEd(sacI) = (endTTs - startTTs)/1000;

        end
    end
end

sacTABLE_Eye0 = table(sac_eyeZEROp , sac_eyeZEROd,...
    'VariableNames',{'SAC_Eye0_Points','SAC_Eye0_DurSecs'});

sacTABLE_Eye1 = table(sac_eyeONEp , sac_eyeONEd,...
    'VariableNames',{'SAC_Eye1_Points','SAC_Eye1_DurSecs'});


end