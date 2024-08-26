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

function [tsTable, fixTab, saccTab, rawTab] = ExtractEyeInfo_MO(eyeMatfile)


load(eyeMatfile,'edfRawStruct');
edfRAW = edfRawStruct;

% Picture and Timestamp Tables
[tsTable] = createTStable(edfRAW);

% Get fixation data
fixTab = table(transpose(edfRAW.Events.Sfix.time),...
    transpose(edfRAW.Events.Sfix.eye),transpose(edfRAW.Events.Efix.eye),...
    transpose(edfRAW.Events.Efix.start),transpose(edfRAW.Events.Efix.end),...
    transpose(edfRAW.Events.Efix.duration),transpose(edfRAW.Events.Efix.posX),...
    transpose(edfRAW.Events.Efix.posY),transpose(edfRAW.Events.Efix.pupilSize),...
    'VariableNames',{'STime','SEye','EEye','EStart','EEnd','Eduration','PosX',...
    'PoxY','PupilS'});

% Get Saccade data
saccTab = table(transpose(edfRAW.Events.Ssacc.time),...
    transpose(edfRAW.Events.Ssacc.eye),transpose(edfRAW.Events.Ssacc.eye),...
    transpose(edfRAW.Events.Esacc.start),transpose(edfRAW.Events.Esacc.end),...
    transpose(edfRAW.Events.Esacc.duration),transpose(edfRAW.Events.Esacc.posX),...
    transpose(edfRAW.Events.Esacc.posY),transpose(edfRAW.Events.Esacc.posXend),...
    transpose(edfRAW.Events.Esacc.posYend),transpose(edfRAW.Events.Esacc.pvel),...
    transpose(edfRAW.Events.Esacc.value1),transpose(edfRAW.Events.Esacc.value2),...
    'VariableNames',{'STime','SEye','EEye','EStart','EEnd','Eduration','PosX',...
    'PosY','PosXend','PosYend','Vel','Value1','Value2'});

% Raw tab
% time, posX, posY, gxvel, gyvel, pa 
rawTab = table(edfRAW.Samples.time,edfRAW.Samples.posX,edfRAW.Samples.posY,...
    edfRAW.Samples.gxvel, edfRAW.Samples.gyvel,edfRAW.Samples.pa,'VariableNames',...
    {'Time','PosX','PosY','VelX','VelY','PupilS'});
end


function [tsTable] = createTStable(edfR)


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

% % picture table
% % Timestamp, picture name, Trial number, 
% mseColPt = contains(mseCol_Char,'TRIAL');
% messPIC = mseCol(mseColPt);
% triIDcont = split(messPIC, ' ');
% triID = cellfun(@(x) str2double(x), triIDcont(:,:,2), "UniformOutput", true);
% timeStamp2 = tsCol(mseColPt);
% 
% fileLOC1 = triIDcont(:,:,3);
% 
% picLOC = cell(length(fileLOC1),1);
% picName = cell(length(fileLOC1),1);
% for pi = 1:length(fileLOC1)
% 
%     [filepath,name,ext] = fileparts(fileLOC1{pi});
% 
%     picLOC{pi} = filepath;
%     picName{pi} = [name , ext];
% 
% end
% 
% picTable = table(transpose(triID),transpose(timeStamp2),picName,...
%     picLOC,'VariableNames',{'TrialNum','timeStamp','Picture','PicLocation'});
end