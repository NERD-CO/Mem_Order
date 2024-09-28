function [] = EyeAnalysis_DataExtract_Gaze_MO(excelLOC , mainLOC, ptID, saveLOC, stimLOC)

% Location of variant.xlsx
cd(excelLOC)
varTable = readtable('Eye_file_information.xlsx');

% Location of eye tracking folders
cd(mainLOC)
[outFOLDS] = getfiles(mainLOC,1,nan);

% Match specific patient ID to location in outFOLDS
ptIdx = find(strcmp(outFOLDS, ptID));

tempCASEd = [mainLOC , filesep , outFOLDS{ptIdx}];
idTab = varTable(matches(varTable.Subject,outFOLDS{ptIdx}),:);
cd(tempCASEd)

% Available MO sessions
avail_MOsessions = ~matches(idTab.EyeData,'NaN');
all_MOsessions = idTab.TaskID(avail_MOsessions);

procFolder = [mainLOC,filesep , ptID, filesep, 'Eye-tracking\Processed'];
behFolder = [mainLOC,filesep , ptID, filesep, 'Behavioral_Data\Raw'];
behFolderP = [mainLOC,filesep , ptID, filesep, 'Behavioral_Data\Processed'];
cd(procFolder)

% Get data files
matMOdirlist = dir('*.mat');
matMOfnames = {matMOdirlist.name};
matMOfparts = split(matMOfnames,{'_','.'});
if ndims(matMOfparts) == 3
    matMOf_ids = matMOfparts(:,:,3);
else
    matMOf_ids = matMOfparts(3);
end

for vi = 1:length(all_MOsessions)

    % Grab session ID
    mo_sessionID = all_MOsessions{vi};

    % Patient answers - need to change per pt - can add input arg for pt file
    cd(behFolder); % Used to create variable 'outData' with TTL values

    % Convert behavioral text files to mat files
    matLIST1 = dir('*.mat');
    matLIST2 = {matLIST1.name};
    switch mo_sessionID
        case 'en'
            moSessionMatFname = matLIST2{contains(matLIST2,'encoding')};
        case 'sc'
            moSessionMatFname = matLIST2{contains(matLIST2,'scene')};
        case 'ti'
            moSessionMatFname = matLIST2{contains(matLIST2,'timeDiscrim')};
    end

    % Load in relevant Picture VARIANT folder
    switch mo_sessionID
        case 'en'
            imageLOC = [stimLOC , filesep , 'A02_Clips'];
        case 'sc'
            imageLOC = [stimLOC , filesep , 'B01_SceneRecogImg'];
        case 'ti'
            imageLOC = [stimLOC , filesep , 'C01_TimeDisImg'];
    end

    % Load respmat
    respMat = load(moSessionMatFname);

    fieldNames = fieldnames(respMat);
    matchedNames = cell(length(fieldNames),1);
    for ii = 1:length(fieldNames)
        if contains(fieldNames{ii,1},'respMat')
            matchedNames{ii,1} = fieldNames{ii,1};
        end
    end

    matchedNames = matchedNames(~cellfun('isempty',matchedNames));
    respMat = respMat.(matchedNames{1,1});

    % Get eye data
    cd(procFolder)
    tmpEye_matFile = matMOfnames{matches(matMOf_ids,mo_sessionID)};
    [tsTable,vidQtable, ~ , ~ , ~, ~, allRAWgxgy_L] =...
        ExtractEyeInfo_EL_MO(tmpEye_matFile,respMat);

    vidQtable.PicLocation = repmat({imageLOC},height(vidQtable),1);

    vidQtable = getPICinfo_MO(vidQtable,mo_sessionID);
    % tsT_L = getTRIinfo(tsT_L);

    % Trial IDs
    trialID = unique(tsTable.trialID(tsTable.trialID ~= 0 & ~isnan(tsTable.trialID)));

    [leftEYE , rightEYE , TTLinfo] = getEYErawEpoch_GAZE_MO(allRAWgxgy_L ,...
        tsTable , trialID , vidQtable , imageLOC, 500, 500 , respMat , mo_sessionID);

    % Convert respMat respValues
    respMat = convertrespVals(respMat);

    % Create table
    task = extractBefore(extractAfter(moSessionMatFname,'_'),'.');
    outTable = generateTable(respMat,task);

    %% Add to structure
    outInfo.(task).table = outTable;
    outInfo.(task).leftEYE = leftEYE;
    outInfo.(task).rightEYE = rightEYE;
    outInfo.(task).tsTable = tsTable;
    outInfo.(task).TTLinfo = TTLinfo;
    outInfo.(task).vidQtable = vidQtable;
    outInfo.(task).respMat = respMat;


end
% subject
saveFname = ['eyeDataGAZE_',ptID,'.mat'];
cd(saveLOC)
save(saveFname,"outInfo",'-v7.3');


%end

end % END MAIN FUNCTION



function [outfiles] = getfiles(dirIN,stage,ftype)

cd(dirIN)
switch stage
    case 1

        foldeS = dir();
        foldeS2 = {foldeS.name};
        foldeS3 = foldeS2(~ismember(foldeS2,{'.','..'}));
        outfiles = foldeS3;
    case 2

        filES = dir(['*.',ftype]);
        filES2 = {filES.name};
        outfiles = filES2;

end


end




% 
% % Trial info
% function [trialINFOtab] = getTRIinfo(inTable)
% 
% ttlIDt = inTable.TTLid;
% 
% ttltrialID = zeros(length(ttlIDt),1);
% trialcount = 0;
% for ti = 1:length(ttltrialID)
%     tmpT = ttlIDt(ti);
%     if tmpT == 1
%         trialcount = trialcount + 1;
%         ttltrialID(ti) = trialcount;
%     else
%         ttltrialID(ti) = trialcount;
%     end
% 
% end
% 
% % start and end
% ttltrialID(ttlIDt == 55) = 0;
% ttltrialID(ttlIDt == 66) = 0;
% 
% trialINFOtab = inTable;
% trialINFOtab.TrialID = ttltrialID;
% 
% 
% end



% Get TS info
function [tsBlk_OUT] = getTSBlock(startI,endI,rawT)

[~, eyeTTL1_i] = min(abs(double(startI) - rawT.TimeTT));
[~, eyeTTL2_i] = min(abs(double(endI) - rawT.TimeTT));

tsBlk_OUT = rawT(eyeTTL1_i:eyeTTL2_i,:);


end


% Clean up Pos
function [cleanPOS] = cleanUPpos(posIN , type , xbounds , ybounds)

% long values
% posL1 = posT2(posT2(:,1) > 1,:);
% posL2 = posL1(posL1(:,2) > 1,:);

switch type
    case 1
        cleanPOS = double(posIN);
        cleanPOS(cleanPOS > 10000) = nan;
    case 2
        inframeTEST = zeros(height(posIN),2,'logical');
        tmpDOUBLE = double(posIN);
        inframeTEST(:,1) = tmpDOUBLE(:,1) > 0 & tmpDOUBLE(:,1) < 768;
        inframeTEST(:,2) = tmpDOUBLE(:,2) > 0 & tmpDOUBLE(:,2) < 1024;
        cleanPOS = all(inframeTEST,2);
    case 3
        inframeTEST = zeros(height(posIN),2,'logical');
        tmpDOUBLE = double(posIN);
        inframeTEST(:,1) = tmpDOUBLE(:,2) > xbounds(1) &...
            tmpDOUBLE(:,2) < xbounds(2);
        inframeTEST(:,2) = tmpDOUBLE(:,1) > ybounds(1) &...
            tmpDOUBLE(:,1) < ybounds(2);
        cleanPOS = all(inframeTEST,2);

        % plot(cleanPOS(:,2),cleanPOS(:,1));ylim([0 768]);xlim([0 1024])
        % xline(xbounds)
        % yline(ybounds)

end

end



function [outEye1 , outEye2] = createEYEtable(pos1,pos2,posF1,posF2,...
    posP1,posP2,timeIN)

for ei = 1:2
    switch ei
        case 1
            eyedataR = {pos1};
            eyeCenR = {mean(pos1,'omitnan')};
            eyeSDR = {std(pos1,'omitnan')};

            eyedataF = {pos1(posF1,:)};
            eyeCenF = {mean(pos1(posF1,:),'omitnan')};
            eyeSDF = {std(pos1(posF1,:),'omitnan')};

            eyedataP = {pos1(posP1,:)};
            eyeCenP = {mean(pos1(posP1,:),'omitnan')};
            eyeSDP = {std(pos1(posP1,:),'omitnan')};

            timeFr = {timeIN(posF1)};
            timePic = {timeIN(posP1)};

            outEye1 = table(eyedataR,eyeCenR,eyeSDR,...
                eyedataF,eyeCenF,eyeSDF,...
                eyedataP,eyeCenP,eyeSDP,...
                timeFr,timePic,'VariableNames',...
                {'gaze_Raw','gaze_Raw_cen','gaze_Raw_sd',...
                'gaze_Scr','gaze_Scr_cen','gaze_Scr_sd',...
                'gaze_Pic','gaze_Pic_cen','gaze_Pic_sd',...
                'gaze_Scr_time','gaze_Pic_time'});

        case 2

            eyedataR = {pos2};
            eyeCenR = {mean(pos2,'omitnan')};
            eyeSDR = {std(pos2,'omitnan')};

            eyedataF = {pos2(posF2,:)};
            eyeCenF = {mean(pos2(posF2,:),'omitnan')};
            eyeSDF = {std(pos2(posF2,:),'omitnan')};

            eyedataP = {pos2(posP2,:)};
            eyeCenP = {mean(pos2(posP2,:),'omitnan')};
            eyeSDP = {std(pos2(posP2,:),'omitnan')};

            timeFr = {timeIN(posF2)};
            timePic = {timeIN(posP2)};
 
            outEye2 = table(eyedataR,eyeCenR,eyeSDR,...
                eyedataF,eyeCenF,eyeSDF,...
                eyedataP,eyeCenP,eyeSDP,...
                timeFr,timePic,'VariableNames',...
                {'gaze_Raw','gaze_Raw_cen','gaze_Raw_sd',...
                'gaze_Scr','gaze_Scr_cen','gaze_Scr_sd',...
                'gaze_Pic','gaze_Pic_cen','gaze_Pic_sd',...
                'gaze_Scr_time','gaze_Pic_time'});

    end


end
end





function [leftEYE , rightEYE , TTL_sInfo] = getEYErawEpoch(rawTIME ,...
    ttlTABLE , trialsOfInt , picTABLE , picLOCation)

%%%%% NEW RAW will go from -400ms [from start][total 500 ITI] to +200ms [from end]
%%%%% NEW TABLE per trial: EVENT , TRIALnum , NLX_T [future] , EYELink_T
%%% TABLE: TRIAL_Num , EVENT_ID , NLX_TTL , EYElink_TTL , EYElink_Int

% Make N by 2 matrix of fieldname + value type
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

% Make table using fieldnames & value types from above
leftEYE = table('Size',[0,size(bothVarNames,1)],...
    'VariableNames', bothVarNames(:,1),...
    'VariableTypes', bothVarNames(:,2));

rightEYE = table('Size',[0,size(bothVarNames,1)],...
    'VariableNames', bothVarNames(:,1),...
    'VariableTypes', bothVarNames(:,2));

TTL_sInfo = cell(numel(unique(ttlTABLE.TrialID(ttlTABLE.TrialID ~=0))),1);
for tttrialir = 1:length(trialsOfInt)
    tmpOtr = trialsOfInt(tttrialir);

    tmpOtrTAB = ttlTABLE(ismember(ttlTABLE.TrialID, tmpOtr),:);

    % 500ms before first TS
    startTS = tmpOtrTAB.timeStamp(tmpOtrTAB.TTLid == 1) - 500;
    % 500ms after last TS
    endTS   = tmpOtrTAB.timeStamp(tmpOtrTAB.TTLid == 6) + 200;

    tmpTable = table;
    tmpTable.TrialNUM = repmat(tmpOtr,height(tmpOtrTAB) + 2,1);
    tmpTable.EventID = [nan ; tmpOtrTAB.TTLid ; nan];
    tmpTable.EventDesc = {'m500ms' ; 'StimOn'; 'StimOff';'QuestionScr';...
        'Response';'DelayEnd';'p200ms'};
    tmpTable.NLXttl = zeros(height(tmpOtrTAB.TTLid) + 2,1,'int64');
    tmpTable.ELNKttl = [startTS ; tmpOtrTAB.timeStamp ; endTS];

    [tsBlk_OUT] = getTSBlock(startTS,endTS,rawTIME);

    middleEvents = zeros(height(tmpOtrTAB),1,'int32');
    for mi = 1:height(tmpOtrTAB)

        [~,elnkLoc] = min(abs(double(tmpOtrTAB.timeStamp(mi)) - tsBlk_OUT.TimeTT));
        middleEvents(mi) = elnkLoc;

    end

    tmpTable.ELNKint = [1 ; middleEvents ; height(tsBlk_OUT)];

    TTL_sInfo{tttrialir} = tmpTable;

    %%% LOAD in IMAGE
    picID = [picLOCation , filesep , picTABLE.CatID{tttrialir},...
        filesep , picTABLE.Picture{tttrialir}];

    [xPictBnds , yPictBnds] = getPICtureDims(picID);

    % EYE 1
    pos_1eye_raw = [tsBlk_OUT.GX_0 , tsBlk_OUT.GY_0];
    pos_1eyeC_raw = cleanUPpos(pos_1eye_raw , 1);
    % IN SCREEN
    pos_1eye_InFr_ind = cleanUPpos(pos_1eye_raw,2);
    % IN Picture
    pos_1eye_InPic_ind = cleanUPpos(pos_1eye_raw,3,xPictBnds,yPictBnds);

    % EYE 2
    pos_2eye_raw = [tsBlk_OUT.GX_1 , tsBlk_OUT.GY_1];
    pos_2eyeC_raw = cleanUPpos(pos_2eye_raw,1);
    % IN SCREEN
    pos_2eye_InFr_ind = cleanUPpos(pos_2eye_raw,2);
    % IN Picture
    pos_2eye_InPic_ind = cleanUPpos(pos_2eye_raw,3,xPictBnds,yPictBnds);

    % USE INDEX OF IN FRAME TO GET POINTs AND TIMEpoints
    
    [left1 , right2] = createEYEtable(pos_1eyeC_raw,pos_2eyeC_raw,...
        pos_1eye_InFr_ind,pos_2eye_InFr_ind,pos_1eye_InPic_ind,...
        pos_2eye_InPic_ind,tsBlk_OUT.TimeTT);
    % left1.catID = old_catID(oTrial);
    % right2.catID = old_catID(oTrial);

    leftEYE(tttrialir,:) = left1;
    rightEYE(tttrialir,:) = right2;

end



end




function [outTABLE2] = fixANDcombineTables(left_learn,right_learn,left_recog,...
    right_recog, picture_learn, picture_recog, learnTTLsummary,...
    recogTTLsummary, learnTTLinfo, recogTTLinfo , recogInfo)

% Left Learn
leftlearnVars = left_learn.Properties.VariableNames;
leftlearnVarsN = cellfun(@(x) ['Left_L_', x], leftlearnVars, 'UniformOutput',false);
left_learn = renamevars(left_learn, leftlearnVars, leftlearnVarsN);

% Left Recog
leftRecogVars = left_recog.Properties.VariableNames;
leftRecoVarsN = cellfun(@(x) ['Left_R_', x], leftRecogVars, 'UniformOutput',false);
left_recog = renamevars(left_recog, leftRecogVars, leftRecoVarsN);

% Right Learn
rightlearnVars = right_learn.Properties.VariableNames;
rightlearnVarsN = cellfun(@(x) ['Right_L_', x], rightlearnVars, 'UniformOutput',false);
right_learn = renamevars(right_learn, rightlearnVars, rightlearnVarsN);

% Right Recog
rightRecogVars = right_recog.Properties.VariableNames;
rightRecogVarsN = cellfun(@(x) ['Right_R_', x], rightRecogVars, 'UniformOutput',false);
right_recog = renamevars(right_recog, rightRecogVars, rightRecogVarsN);

% Picture learn
learnPicVars = picture_learn.Properties.VariableNames;
learnPicVarsN = cellfun(@(x) ['Learn_', x], learnPicVars, 'UniformOutput',false);
picture_learn = renamevars(picture_learn, learnPicVars, learnPicVarsN);

% Picture recog
recogPicVars = picture_recog.Properties.VariableNames;
recogPicVarsN = cellfun(@(x) ['Recog_', x], recogPicVars, 'UniformOutput',false);
picture_recog = renamevars(picture_recog, recogPicVars, recogPicVarsN);

outTABLE1 = [left_learn , left_recog , right_learn , right_recog ,...
    picture_learn , picture_recog];

outTABLE1.LearnTTL = learnTTLsummary;
outTABLE1.RecogTTL = recogTTLsummary;
outTABLE1.LearnTTLplus = learnTTLinfo;
outTABLE1.RecogTTLplus = recogTTLinfo;

outTABLE2 = [outTABLE1 , recogInfo];


end




function [outTABLE] = convertRAW2trial(inTTLtab)


ttlINFO = inTTLtab.taskinformation;

indexLOG = matches(ttlINFO.TTLvalue,'1');
indexLOC = find(indexLOG);
% trialN = sum(indexLOG);
% trialVec = 1:trialN;
trialCell = cell(sum(indexLOG),1);
trialInVec = zeros(sum(~matches(ttlINFO.TTLvalue,'66')),1);

for ii = 1:numel(indexLOC)

    if ii ~= length(indexLOC)
        tmpStart = indexLOC(ii);
        tmpEnd = indexLOC(ii + 1) - 1;
        trialInVec(tmpStart:tmpEnd) = ii;
        tmpTAB = ttlINFO(tmpStart:tmpEnd,:);
        tmpTAB.TrialID = trialInVec(tmpStart:tmpEnd);
        trialCell{ii} = tmpTAB;
    else
        tmpStart = indexLOC(ii);
        tmpEnd = sum(~matches(ttlINFO.TTLvalue,'66'));
        trialInVec(tmpStart:tmpEnd) = ii;
        tmpTAB = ttlINFO(tmpStart:tmpEnd,:);
        tmpTAB.TrialID = trialInVec(tmpStart:tmpEnd);
        trialCell{ii} = tmpTAB;
    end
end

outTABLE = trialCell;


end




function [xBounds , yBounds] = getPICtureDims(imageLOCATION)

img = imread(imageLOCATION);
[xRow,yCol] = size(img,[1,2]);
screenWidth = 1024;
screenHeight = 768;
% Image dimensions
imageWidth = xRow;
imageHeight = yCol;

% Calculate centered coordinates -- top left corner of centered image
x = (screenWidth - imageWidth) / 2;
y = (screenHeight - imageHeight) / 2;

xBounds = round([x x+imageWidth]);
yBounds = round([y y+imageHeight]);

end