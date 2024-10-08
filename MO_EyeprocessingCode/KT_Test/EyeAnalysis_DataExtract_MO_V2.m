function [] = EyeAnalysis_DataExtract_MO_V2(excelLOC , mainLOC, ptID, saveLOC)

%% Location of variant.xlsx
cd(excelLOC)
varTable = readtable('Eye_file_information.xlsx');

%% Location of eye tracking folders
cd(mainLOC)
[outFOLDS] = getfiles(mainLOC,1,nan);

%% Match specific patient ID to location in outFOLDS
ptIdx = find(strcmp(outFOLDS, ptID));

tempCASEd = [mainLOC , filesep , outFOLDS{ptIdx}];
idTab = varTable(matches(varTable.Subject,outFOLDS{ptIdx}),:);
cd(tempCASEd)

%% Available MO sessions
avail_MOsessions = ~matches(idTab.EyeData,'NaN');
all_MOsessions = idTab.TaskID(avail_MOsessions);

%% Folders
procFolder = [mainLOC,filesep , ptID, filesep, 'Eye-tracking\Processed'];
behFolder = [mainLOC,filesep , ptID, filesep, 'Behavioral_Data\Raw'];
behFolderP = [mainLOC,filesep , ptID, filesep, 'Behavioral_Data\Processed'];
cd(procFolder)

%% Get data files
matMOdirlist = dir('*.mat');
matMOfnames = {matMOdirlist.name};
matMOfparts = split(matMOfnames,{'_','.'});
if ndims(matMOfparts) == 3
    matMOf_ids = matMOfparts(:,:,3);
else
    matMOf_ids = matMOfparts(3);
end

%% Create structure
moSessionSall = struct;

%% Loop through sessions
for vi = 1:length(all_MOsessions) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Grab session ID
    mo_sessionID = all_MOsessions{vi};
    % Patient answers - need to change per pt - can add input arg for pt file
    cd(behFolder); % Used to create variable 'outData' with TTL values

    %% Convert behavioral text files to mat files
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
        case 'ti'
            moSessionTxtFname = txtLIST2{contains(txtLIST2,'timeDiscrim')};
            moSessionMatFname = matLIST2{contains(matLIST2,'timeDiscrim')};
    end

    %% Load respMat
    %load(moSessionMatFname,'respMat');
    respMat = load(moSessionMatFname);

    %% Find field containing respMat data
    fieldNames = fieldnames(respMat);
    matchedNames = cell(length(fieldNames),1);
    for ii = 1:length(fieldNames)
        if contains(fieldNames{ii,1},'respMat')
            matchedNames{ii,1} = fieldNames{ii,1};
        end
    end

    %% Remove empty cells
    matchedNames = matchedNames(~cellfun('isempty',matchedNames));

    %% Grab respMat
    respMat = respMat.(matchedNames{1,1});

    %% Get eye data
    cd(procFolder)
    tmpEye_matFile = matMOfnames{matches(matMOf_ids,mo_sessionID)};
    [tsTable,vidQtable,~,~,sacTABLE_Eye0,sacTABLE_Eye1,allgazeTAB] =...
        ExtractEyeInfo_EL_MO(tmpEye_matFile,respMat);

    %% Trial IDs
    trialID = unique(tsTable.trialID(tsTable.trialID ~= 0 & ~isnan(tsTable.trialID)));

    %% Create eye epochs
    % JAT added respmat input ----- 9/14/2024
    [leftEYE , rightEYE , TTLinfo] = getEYErawEpoch(allgazeTAB , tsTable , trialID, 500, 500 , respMat); %%%%%%%%%%%%%%%%%%%

    %% Add NLX info
    [TTLinfo] = modifyTTLinfo(TTLinfo,ptID,moSessionMatFname,respMat); %%%%%%%%%%%%%%%

    %% Patient ID
    patientID = idTab.Subject{1};

    %% Read function and load in new .mat file
    cd(behFolder);
    % [behavFILE_learn] = MO_TxttoMat(moSessionTxtFname,patientID,mo_sessionID, behFolderP);
    % load(behavFILE_learn, 'outData');
    % tmpMOTTL = outData;
    % cd(mainLOC);
    % tmpMOTTLsummary = convertRAW2trial(tmpMOTTL);

    %% Convert respMat respValues
    respMat = convertrespVals(respMat);

    %% Create table
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
    outInfo.(task).sacTABLE_Eye0 = sacTABLE_Eye0;
    outInfo.(task).sacTABLE_Eye1 = sacTABLE_Eye1;

    disp(['Finished SESSION # ',num2str(vi)])

end

saveFname = ['eyeData_',ptID,'.mat'];
cd(saveLOC)
save(saveFname,"outInfo",'-v7.3');

end