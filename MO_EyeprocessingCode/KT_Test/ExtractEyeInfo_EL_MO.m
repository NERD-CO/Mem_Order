function [tsTable,vidQuestable,fixTABLE_Eye0,fixTABLE_Eye1,sacTABLE_Eye0,sacTABLE_Eye1,allgazePosPupsTAB] = ExtractEyeInfo_EL_MO(eyeMatfile , behavOUTPUT)

%% function [tsTable,vidQuestable,fixTABLE_Eye0,fixTABLE_Eye1,sacTABLE_Eye0,sacTABLE_Eye1,allgazePosPupsTAB] = ExtractEyeInfo_EL_MO(eyeMatfile , behavOUTPUT)
% Extract eyetracking data from converted .edf to .mat file
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

%% Load the data
load(eyeMatfile,'edfRawStruct');
edfRAW = edfRawStruct;

%% Picture and timestamp tables
[tsTable, vidQuestable] = createTStable(edfRAW , behavOUTPUT);

%%

%% Grab data
edfRAWEvent = struct2table(edfRAW.FEVENT);
gxALL = edfRAW.FSAMPLE.gx;
gyALL = edfRAW.FSAMPLE.gy;
timeALL = edfRAW.FSAMPLE.time;
posX = edfRAW.FSAMPLE.px;
posY = edfRAW.FSAMPLE.py;
velX = edfRAW.FSAMPLE.gxvel;
velY = edfRAW.FSAMPLE.gyvel;
pupilS = edfRAW.FSAMPLE.pa;

%% Convert to table
allgazePosPupsTAB = array2table([transpose(timeALL) , transpose(posX) , transpose(posY) ,...
    transpose(velX), transpose(velY) , transpose(gxALL),...
    transpose(gyALL) , transpose(pupilS)],'VariableNames',...
    {'Time','PosX_0','PosX_1','PosY_0','PosY_1','VelX_0','VelX_1','VelY_0',...
    'VelY_1','GX_0','GX_1','GY_0','GY_1','PupilS_0','PupilS_1'});

%% Set eye IDs
eyeIDs = [0 1];
allEYEs = edfRAWEvent.eye;
eyesOUT = struct;

%% Loop over eyes for fixes
for eyeI = 1:2
    tmpEYE = eyeIDs(eyeI);
    eyeTAB = edfRAWEvent(allEYEs == tmpEYE,:);

    % Get FIX start %
    endFIX = matches(eyeTAB.codestring,'ENDFIX');
    eFIXtab = eyeTAB(endFIX,:);

    eyesOUT.(['eye_',num2str(tmpEYE)]).Fixes = [eFIXtab.sttime,...
        eFIXtab.entime];
end

%% Label start and end times
eyesOUT.eye_0.Fixes = array2table(eyesOUT.eye_0.Fixes,'VariableNames',...
    {'StartTime','EndTime'});
eyesOUT.eye_1.Fixes = array2table(eyesOUT.eye_1.Fixes,'VariableNames',...
    {'StartTime','EndTime'});

%% Loop over eyes to add saccades
for eyeI = 1:2

    tmpEYE = eyeIDs(eyeI);
    eyeTAB = edfRAWEvent(allEYEs == tmpEYE,:);

    % Get FIX start
    endSAC = matches(eyeTAB.codestring,'ENDSACC');
    eSACtab = eyeTAB(endSAC,:);

    eyesOUT.(['eye_',num2str(tmpEYE)]).Saccades = [eSACtab.sttime,...
        eSACtab.entime];
end

%% Label start and end times
eyesOUT.eye_0.Saccades = array2table(eyesOUT.eye_0.Saccades,'VariableNames',...
    {'StartTime','EndTime'});
eyesOUT.eye_1.Saccades = array2table(eyesOUT.eye_1.Saccades,'VariableNames',...
    {'StartTime','EndTime'});

%% Process raw gx and gy positions
[fixTABLE_Eye0,fixTABLE_Eye1,sacTABLE_Eye0,sacTABLE_Eye1] = getFIXSACdata(eyesOUT,gxALL,gyALL,timeALL);

end