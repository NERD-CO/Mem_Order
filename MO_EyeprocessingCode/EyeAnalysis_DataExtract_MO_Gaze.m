function [] = EyeAnalysis_DataExtract_MO_Gaze(excelLOC , mainLOC, ptID, saveLOC)
% 'The events coorespond to the TTL markers for each trial. ', ...
%     'For the learning trials, the TTL markers are the following: 55 = start of the experiment, ', ...
%     '1 = stimulus ON, 2 = stimulus OFF, 3 = Question Screen Onset, ', ...
%     '20 = Yes (21 = NO) during learning, 6 = End of Delay after Response, ', ...
%     '66 = End of Experiment. For the recognition trials, the TTL markers are the following: ', ...
%     '55 = start of experiment, 1 = stimulus ON, 2 = stimulus OFF, 3 = Question Screen Onset, ' ...
%     '31:36 = Confidence (Yes vs. No) response, 66 = End of Experiment'

% Loop through files - determine if learn or retrieve

% learning = [55, 1, 2, 3 20, 21, 6, 66];
% recog = [55, 1, 2, 3, 31:36, 66];

% Location of variant.xlsx
cd(excelLOC)
varTable = readtable('variantLIST.xlsx');
% Location of eye tracking folders
cd(mainLOC)
[outFOLDS] = getfiles(mainLOC,1,nan);

% Match specific patient ID to location in outFOLDS
ptIdx = find(strcmp(outFOLDS, ptID));

tempCASEd = [mainLOC , filesep , outFOLDS{ptIdx}];
idTab = varTable(matches(varTable.Subject,outFOLDS{ptIdx}),:);
cd(tempCASEd)

allvars = unique(idTab.Variant);
variantS = struct;
for vi = 1:length(allvars)

    % Get Variant and Condition

    variantTAB = idTab(ismember(idTab.Variant,allvars(vi)),:);
    learn_matFile1 = variantTAB.MatFile{matches(variantTAB.Block,'learn')};
    learn_matFile = replace(learn_matFile1,'eyetrack','eyetrackELF');

    % Load in relevant Picture VARIANT folder
    switch variantTAB.Variant(1)
        case 1
            imageLOC = 'D:\Dropbox\LearningRecognitionECoG\Stimuli\newolddelay';
        case 2
            imageLOC = 'D:\Dropbox\LearningRecognitionECoG\Stimuli\newolddelay2';
        case 3
            imageLOC = 'D:\Dropbox\LearningRecognitionECoG\Stimuli\newolddelay3';
    end


    [tsT_L, picT_L, eye0raw_L , ~ , eye1raw_L, ~, allRAWgxgy_L] = ExtractEyeInfo_EL(learn_matFile);

    picT_L = getPICinfo(picT_L);
    tsT_L = getTRIinfo(tsT_L);

    trialID_L = unique(tsT_L.TrialID(tsT_L.TrialID ~= 0));

    [leftEYE_Learn , rightEYE_Learn , learnTTLinfo] = getEYErawEpoch(allRAWgxgy_L ,...
        tsT_L , trialID_L , picT_L , imageLOC);

    % Patient answers - need to change per pt - can add input arg for pt file
    cd(tempCASEd); % Used to create variable 'outData' with TTL values

    % Convert behavioral txt files to mat files
    learnTXT = variantTAB.Behavior{matches(variantTAB.Block,'learn')};

    % Recog file
    block = 'learn'; %Will be either 'learn' or 'recog'
    patientID = idTab.Subject{1};

    % Run function and load in new .mat file
    [behavFILE_learn] = NewOldTxttoMat_v2(learnTXT,patientID,vi,block, tempCASEd);
    load(behavFILE_learn, 'outData') %loc and file name as input
    learnTTL = outData;
    cd(mainLOC);
    learnTTLsummary = convertRAW2trial(learnTTL);

    ttlValuesLearn = learnTTL.taskinformation.TTLvalue;
    yesNoLearn = ttlValuesLearn(matches(ttlValuesLearn,{'20','21'}));


    % RECOGNITION PROCESSING ----------------------------------------------

    % load in information for image paths and whether new/old in recog
    cd(excelLOC)
    load('newOld_stimID_all.mat','stimAll');

    % add in compareAns code here
    % Extract variables for ground truth old vs new- 1=present in learn block
    groundTruthRecog_var1 = stimAll.stimNewOld_var1;
    groundTruthRecog_var2 = stimAll.stimNewOld_var2;
    groundTruthRecog_var3 = stimAll.stimNewOld_var3;

    % Patient answers - need to change per pt - can add input arg for pt file
    cd(tempCASEd); % Used to create variable 'outData' with TTL values

    % Convert behavioral txt files to mat files
    recogTXT = variantTAB.Behavior{matches(variantTAB.Block,'recog')};

    % Recog file
    block = 'recog'; %Will be either 'learn' or 'recog'
    patientID = idTab.Subject{1};

    % Run function and load in new .mat file
    [behavFILE_recog] = NewOldTxttoMat_v2(recogTXT,patientID,vi,block, tempCASEd);
    load(behavFILE_recog, 'outData') %loc and file name as input
    recogTTL = outData;
    cd(mainLOC);
    recogTTLsummary = convertRAW2trial(recogTTL);

    % Extract answers from outData - 31:36 = Confidence (Yes vs. No) response
    % 31,32,33=no,new & 34,35,36=yes,old
    ttlValues = str2double(recogTTL.taskinformation.TTLvalue);
    confRatings = ttlValues(ttlValues(:,1) >= 31 & ttlValues(:,1) <= 36,:);

    % Logical array of confRatings: GREATER than 34
    confRatings_logical = logical(confRatings >= 34); %1=yes,old & 0=no,new

    confVal_1 = logical(confRatings == 31); %1=31, 0= not 31
    confVal_2 = logical(confRatings == 32); %1=32, 0= not 32
    confVal_3 = logical(confRatings == 33); %1=33, 0= not 33
    confVal_4 = logical(confRatings == 34); %1=34, 0= not 34
    confVal_5 = logical(confRatings == 35); %1=35, 0= not 35
    confVal_6 = logical(confRatings == 36); %1=36, 0= not 36

    confVal_34 = logical(confRatings == 33 | confRatings == 34); %1=33-34, 0=not 33-34
    confVal_12 = logical(confRatings == 31 | confRatings == 32); %1=31-32, 0=not 31-32
    confVal_56 = logical(confRatings == 35 | confRatings == 36); %1=35-36, 0=not 35-36
    confVal_1256 = ~logical(confVal_34);                         %1=31+32+35+36, 0=33+34
    confVal_16 = logical(confRatings == 31 | confRatings == 36); %1=31+36, 0=32+33+34+35

    % Can compare ground truth to patient answers
    % 1=yes,old & 0=no,new
    % Combine arrays into a table
    switch allvars(vi)
        case 1
            groundTruthRecog = groundTruthRecog_var1;
        case 2
            groundTruthRecog = groundTruthRecog_var2;
        case 3
            groundTruthRecog = groundTruthRecog_var3;
    end

    % CONSIDER RE-NAMING RECOG INFO
    recogINFO = table(confRatings, yesNoLearn , confRatings_logical, groundTruthRecog, ...
        confVal_1, confVal_2, confVal_3, confVal_4, confVal_5, confVal_6,...
        confVal_34, confVal_12, confVal_56, confVal_1256, confVal_16, 'VariableNames', ...
        {'RecogResp','LearnResp','confRatings', 'groundTruth', 'New_sure', 'New_less_sure', 'New_unsure',...
        'Old_unsure', 'Old_less_sure', 'Old_sure', 'All_unsure', ...
        'New_sureLess_sure', 'Old_sureLess sure', 'All_sureLess_sure', 'All_sure'});

    cd(tempCASEd)
    recog_matFile1 = variantTAB.MatFile{matches(variantTAB.Block,'recog')};
    recog_matFile = replace(recog_matFile1,'eyetrack','eyetrackELF');

    [tsT_R, picT_R, eye0raw_R , ~ , eye1raw_R, ~, allRAWgxgy_R] = ExtractEyeInfo_EL(recog_matFile);

    picT_R = getPICinfo(picT_R);
    tsT_R = getTRIinfo(tsT_R);

    trialID_R = unique(tsT_R.TrialID(tsT_R.TrialID ~= 0));

    [leftEYE_Recog , rightEYE_Recog , recogTTLinfo] = getEYErawEpoch(allRAWgxgy_R ,...
        tsT_R , trialID_R , picT_R , imageLOC);

    % ADD SUMMARY TTLs , ADD learnYES_NO
    [outTABLE] = fixANDcombineTables(leftEYE_Learn,rightEYE_Learn,leftEYE_Recog,...
        rightEYE_Recog , picT_L, picT_R, learnTTLsummary, recogTTLsummary,...
        learnTTLinfo, recogTTLinfo , recogINFO);

    varianTNum = ['var',num2str(allvars(vi))];

    variantS.(varianTNum).dataTable = outTABLE;

end
% subject
saveFname = ['eyeDataGAZE_',outFOLDS{ptIdx},'.mat'];
cd(saveLOC)
save(saveFname,"variantS");


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


% Pic info
function [picINFOtab] = getPICinfo(inTable)

picLOCa = inTable.PicLocation;
fparTs = split(picLOCa,'\');
catID = fparTs(:,11);
% catSubn = cellfun(@(x) str2double(x(1)), numCAT, 'UniformOutput',true);
% catID = cellfun(@(x) x(2:end), numCAT, 'UniformOutput',false);
picJPG = inTable.Picture;
picNUMs = cellfun(@(x) split(x,'.'), picJPG, 'UniformOutput',false);
picNUM = cellfun(@(x) str2double(x{1}), picNUMs, 'UniformOutput',true);
picINFOtab = inTable;
picINFOtab.CatID = catID;
picINFOtab.PicNUM = picNUM;

% tmpCombine = cell(length(catSubn),1);
% for pi = 1:length(catSubn)
%     tmpCombine{pi} = [num2str(catSubn(pi)),'.',num2str(picNUM(pi))];
% end
% picINFOtab.CatPICid = tmpCombine;

end


% Trial info
function [trialINFOtab] = getTRIinfo(inTable)

ttlIDt = inTable.TTLid;

ttltrialID = zeros(length(ttlIDt),1);
trialcount = 0;
for ti = 1:length(ttltrialID)
    tmpT = ttlIDt(ti);
    if tmpT == 1
        trialcount = trialcount + 1;
        ttltrialID(ti) = trialcount;
    else
        ttltrialID(ti) = trialcount;
    end

end

% start and end
ttltrialID(ttlIDt == 55) = 0;
ttltrialID(ttlIDt == 66) = 0;

trialINFOtab = inTable;
trialINFOtab.TrialID = ttltrialID;


end



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