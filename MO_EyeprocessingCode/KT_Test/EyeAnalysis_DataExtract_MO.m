function [] = EyeAnalysis_DataExtract_MO(excelLOC , mainLOC, ptID, saveLOC)

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

avail_MOsessions = ~matches(idTab.EyeData,'NaN');
all_MOsessions = idTab.TaskID(avail_MOsessions);
% all_MOfiles = idTab.EyeData(avail_MOsessions);

procFolder = [mainLOC,filesep , ptID, filesep, 'Eye-tracking\Processed'];
behFolder = [mainLOC,filesep , ptID, filesep, 'Behavioral_Data\Raw'];
behFolderP = [mainLOC,filesep , ptID, filesep, 'Behavioral_Data\Processed'];
cd(procFolder)

matMOdirlist = dir('*.mat');
matMOfnames = {matMOdirlist.name};

matMOfparts = split(matMOfnames,{'_','.'});

matMOf_ids = matMOfparts(:,:,3);

moSessionSall = struct;
for vi = 1:length(all_MOsessions)

    mo_sessionID = all_MOsessions{vi};
    % Patient answers - need to change per pt - can add input arg for pt file
    cd(behFolder); % Used to create variable 'outData' with TTL values

    % Convert behavioral txt files to mat files
    txtLIST1 = dir('*.txt');
    txtLIST2 = {txtLIST1.name};
    matLIST1 = dir('*.mat');
    matLIST2 = {matLIST1.name};
    switch mo_sessionID
        case 'en'
            moSessionTxtFname = txtLIST2{contains(txtLIST2,'encoding')};
            moSessionMatFname = matLIST2{contains(matLIST2,'encoding')};
        case 'sc'
            moSessionTxtFname = txtLIST2{contains(txtLIST2,'scene')};
            moSessionMatFname = matLIST2{contains(matLIST2,'scene')};
        case 'td'
            moSessionTxtFname = txtLIST2{contains(txtLIST2,'timeDiscrim')};
            moSessionMatFname = matLIST2{contains(matLIST2,'timeDiscrim')};
    end

    load(moSessionMatFname,'respMat');

    cd(procFolder)
    % Get Condition
    tmpEye_matFile = matMOfnames{matches(matMOf_ids,mo_sessionID)};
    [tsTable, vidQtable,fixTABLE_Eye0 ,...
        fixTABLE_Eye1 ,...
        sacTABLE_Eye0 ,...
        sacTABLE_Eye1 ,...
        allgazeTAB]  = ExtractEyeInfo_EL_MO(tmpEye_matFile,respMat); 

    trialID = unique(tsTable.TrialID(tsTable.TrialID ~= 0));

    [leftEYE , rightEYE , TTLinfo] = getEYErawEpoch(allgazeTAB , tsTable , trialID); %%%%%%%%%%%%%%%%%%%%% creates tmpTable
    %[leftEYE , rightEYE , TTLinfo] = getEYErawEpoch_V2(allgazeTAB , tsTable , trialID, 400, 200);

    patientID = idTab.Subject{1};

    % Run function and load in new .mat file
    cd(behFolder); 
    [behavFILE_learn] = MO_TxttoMat(moSessionTxtFname,patientID,mo_sessionID, behFolderP); 
    load(behavFILE_learn, 'outData'); %loc and file name as input
    tmpMOTTL = outData;
    cd(mainLOC);
    tmpMOTTLsummary = convertRAW2trial(tmpMOTTL); 


    % CONSIDER RE-NAMING RECOG INFO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% TALK TO KEVIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    recogINFO = table(confRatings, yesNoLearn , confRatings_logical, groundTruthRecog, ...
        confVal_1, confVal_2, confVal_3, confVal_4, confVal_5, confVal_6,...
        confVal_34, confVal_12, confVal_56, confVal_1256, confVal_16, 'VariableNames', ...
        {'RecogResp','LearnResp','confRatings', 'groundTruth', 'New_sure', 'New_less_sure', 'New_unsure',...
        'Old_unsure', 'Old_less_sure', 'Old_sure', 'All_unsure', ...
        'New_sureLess_sure', 'Old_sureLess sure', 'All_sureLess_sure', 'All_sure'});

    % ADD SUMMARY TTLs , ADD learnYES_NO
    [outTABLE] = fixANDcombineTables(leftEYE_Learn,rightEYE_Learn,leftEYE_Recog,...
        rightEYE_Recog , picT_L, picT_R, tmpMOTTLsummary, recogTTLsummary,...
        learnTTLinfo, recogTTLinfo , recogINFO);

    varianTNum = ['var',num2str(allvars(vi))];

    variantS.(varianTNum).dataTable = outTABLE;

end
% subject
saveFname = ['eyeData_',outFOLDS{ptIdx},'.mat'];
cd(saveLOC)
save(saveFname,"variantS");

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

[~, eyeTTL1_i] = min(abs(double(startI) - rawT.Time));
[~, eyeTTL2_i] = min(abs(double(endI) - rawT.Time));

tsBlk_OUT = rawT(eyeTTL1_i:eyeTTL2_i,:);


end


% Clean up Pos
function [cleanPOS] = cleanUPpos(posIN)

posS1 = smoothdata(posIN(:,1),'gaussian',40);
posS2 = smoothdata(posIN(:,2),'gaussian',40);

cleanPOS = [posS1 , posS2];

end



function [outEye1 , outEye2] = createEYEtable(ps1,ps2,pos1,pos2)

for ei = 1:2
    switch ei
        case 1
            eyedata = {pos1};
            eyeCen = {mean(pos1)};
            eyeSD = {std(pos1)};
            Q3pos = quantile(pos1,0.75);
            Q1pos = quantile(pos1,0.25);
            eyeCD = (Q3pos - Q1pos) / (Q3pos + Q1pos);
            eyedist = {pdist2(eyeCen{1},eyedata{1},'euclidean')};
            pupdata = {ps1};
            pupCen = mean(ps1);
            pupSD = std(ps1);
            Q3pup = quantile(ps1,0.75);
            Q1pup = quantile(ps1,0.25);
            pupCD = (Q3pup - Q1pup) / (Q3pup + Q1pup);

            outEye1 = table(eyedata,eyeCen,eyeSD,eyeCD,eyedist,...
                pupdata,pupCen,pupSD,pupCD,'VariableNames',...
                {'oT_posit_raw','oT_posit_cen','oT_posit_sd','oT_posit_cd',...
                'oT_posit_dist','oT_pupilS_raw','oT_pupilS_mean','oT_pupilS_sd',...
                'oT_pupilS_cd'});


        case 2

            eyedata = {pos2};
            eyeCen = {mean(pos2)};
            eyeSD = {std(pos2)};
            Q3pos = quantile(pos2,0.75);
            Q1pos = quantile(pos2,0.25);
            eyeCD = (Q3pos - Q1pos) / (Q3pos + Q1pos);
            eyedist = {pdist2(eyeCen{1},eyedata{1},'euclidean')};
            pupdata = {ps2};
            pupCen = mean(ps2);
            pupSD = std(ps2);
            Q3pup = quantile(ps2,0.75);
            Q1pup = quantile(ps2,0.25);
            pupCD = (Q3pup - Q1pup) / (Q3pup + Q1pup);

            outEye2 = table(eyedata,eyeCen,eyeSD,eyeCD,eyedist,...
                pupdata,pupCen,pupSD,pupCD,'VariableNames',...
                {'oT_posit_raw','oT_posit_cen','oT_posit_sd','oT_posit_cd',...
                'oT_posit_dist','oT_pupilS_raw','oT_pupilS_mean','oT_pupilS_sd',...
                'oT_pupilS_cd'});

    end


end
end





function [leftEYE , rightEYE , TTL_sInfo] = getEYErawEpoch(rawTIME , ttlTABLE , trialsOfInt)

%%%%% NEW RAW will go from -400ms [from start][total 500 ITI] to +200ms [from end]
%%%%% NEW TABLE per trial: EVENT , TRIALnum , NLX_T [future] , EYELink_T
%%% TABLE: TRIAL_Num , EVENT_ID , NLX_TTL , EYElink_TTL , EYElink_Int

% Make N by 2 matrix of fieldname + value type
bothVarNames = [['oT_posit_raw', "cell"]; ...
    ['oT_posit_cen', "cell"]; ...
    ['oT_posit_sd', "cell"]; ...
    ['oT_posit_cd', "double"]; ...
    ['oT_posit_dist', "cell"]; ...
    ['oT_pupilS_raw', "cell"];...
    ['oT_pupilS_mean', "double"];...
    ['oT_pupilS_sd', "double"];...
    ['oT_pupilS_cd', "double"]];

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
    startTS = tmpOtrTAB.timeStamp(tmpOtrTAB.TTLid == 1) - 400;
    % 500ms after last TS
    endTS   = tmpOtrTAB.timeStamp(tmpOtrTAB.TTLid == 3) + 2000;

    tmpTable = table;
    tmpTable.TrialNUM = repmat(tmpOtr,height(tmpOtrTAB) + 2,1);
    tmpTable.EventID = [nan ; tmpOtrTAB.TTLid ; nan];

    if height(tmpTable.EventID) == 6

        tmpTable.EventDesc = {'m400ms' ; 'FixCross'; 'ClipOn';'ClipOff';...
            'Probe';'p2000ms'};

    else

        keyboard; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end
    tmpTable.NLXttl = zeros(height(tmpOtrTAB.TTLid) + 2,1,'int64');
    tmpTable.ELNKttl = [startTS ; tmpOtrTAB.timeStamp ; endTS];

    [tsBlk_OUT] = getTSBlock(startTS,endTS,rawTIME);

    middleEvents = zeros(height(tmpOtrTAB),1,'int32');
    for mi = 1:height(tmpOtrTAB)

        [~,elnkLoc] = min(abs(double(tmpOtrTAB.timeStamp(mi)) - tsBlk_OUT.Time));
        middleEvents(mi) = elnkLoc;

    end

    tmpTable.ELNKint = [1 ; middleEvents ; height(tsBlk_OUT)]; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    TTL_sInfo{tttrialir} = tmpTable;

    pupilS_1 = tsBlk_OUT.PupilS_0(:,1);
    pupilS_2 = tsBlk_OUT.PupilS_1(:,1);

    % NEED TO clean up PupilSize

    % Clean up - nans and high values
    pos_1 = [tsBlk_OUT.GX_0(:,1) , tsBlk_OUT.GY_0(:,1)];
    pos_1c = cleanUPpos(pos_1);
    pos_2 = [tsBlk_OUT.GX_1(:,1) , tsBlk_OUT.GY_1(:,1)];
    pos_2c = cleanUPpos(pos_2);

    [left1 , right2] = createEYEtable(double(pupilS_1),double(pupilS_2),pos_1c,pos_2c);
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

indexLOG = matches(ttlINFO.TTLvalue,'11'); % change from 1 to 11 - KT 8/27/24
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
