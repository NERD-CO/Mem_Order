
% STEP 1: Change path directories. Choose the NAS case but update the paths
% the first time you use it to make sure they are correct.

getPCname = getenv('COMPUTERNAME');

switch getPCname
    case 'DESKTOP-FAGRV5G'
   
    case 'DESKTOP-I5CPDO7'
        % Code location
        codeLocation = 'C:\Users\Admin\Documents\Github\Mem_Order'; %%%%%%%%%%%%%%%%%
        addpath(codeLocation)
        % Mex file for converting eye tracking data
        edfMexELloc = 'C:\Users\Admin\Documents\Github\Mem_Order\MO_EyeprocessingCode\edfmex';
        edfCheck2 = which('edfmex.mexa64');
        if isempty(edfCheck2)
            addpath(genpath(edfMexELloc));
        end

        % OPTIONAL - color maps
        cbrewLOC = 'C:\Users\Admin\Documents\Github\Mem_Order\MO_EyeprocessingCode\cbrewerALL';
        addpath(genpath(cbrewLOC))

        % File with subject information 
        excelLocation = 'Z:\Tyner_K_Projects\MEM_ORDER';
        % Data folder location
        dataLocation = 'Z:\Tyner_K_Projects\MEM_ORDER\DataFolder';
        stimuliLOC = 'Z:\Tyner_K_Projects\MEM_ORDER\StimuliFiles';
        % Save locations for Eye data
        % savePreProcLocation = [dataLocation , filesep , 'eyeDATA'];
        % saveCleanLocation = [savePreProcLocation , filesep , 'cleaned_eyeDATA'];
        cd(codeLocation)

    case 'NERDCO-POSTDOC' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Code location
        codeLocation = 'C:\Users\Kevin_Tyner\Documents\GitHub\Mem_Order\MO_EyeprocessingCode'; %%%%%%%%%%%%%%%%%
        addpath(codeLocation)
        % Mex file for converting eye tracking data
        edfMexELloc = 'C:\Users\Kevin_Tyner\Documents\GitHub\Mem_Order\MO_EyeprocessingCode\edfmex';
        edfCheck2 = which('edfmex.mexa64');
        if isempty(edfCheck2)
            addpath(genpath(edfMexELloc));
        end

        % OPTIONAL - color maps
        cbrewLOC = 'C:\Users\Kevin_Tyner\Documents\GitHub\Mem_Order\MO_EyeprocessingCode\cbrewerALL';
        addpath(genpath(cbrewLOC))

        % File with subject information 
        excelLocation = 'Z:\Tyner_K_Projects\MEM_ORDER';
        % Data folder location
        dataLocation = 'Z:\Tyner_K_Projects\MEM_ORDER\DataFolder';
        stimuliLOC = 'Z:\Tyner_K_Projects\MEM_ORDER\StimuliFiles';
        addpath(genpath('Z:\Tyner_K_Projects\MEM_ORDER\StimuliFiles'))
        % Save locations for Eye data
        % savePreProcLocation = [dataLocation , filesep , 'eyeDATA'];
        % saveCleanLocation = [savePreProcLocation , filesep , 'cleaned_eyeDATA'];
        cd(codeLocation)
 

end

%% STEP 2: Change ptID to be specific to pt 
ptID = 'MW26'; %MW26

%% STEP 3 CONVERT From EDF to MAT
%Extract_Eye_EDF_MO(excelLocation , dataLocation, ptID)
clc

%% STEP 4 KT test Initial_EyeAnalysis function %% originally commented out

savePreProcLocation = [dataLocation , filesep ,ptID , filesep ,...
    'Eye-tracking\Processed\eyeDATA'];

EyeAnalysis_DataExtract_MO_V2(excelLocation, dataLocation, ptID, savePreProcLocation);
clc

%% STEP 4.1 RUN initial eye position function %% Originally commented out

savePreProcLocation = [dataLocation , filesep ,ptID , filesep ,...
    'Eye-tracking\Processed\eyeDATA'];

EyeAnalysis_DataExtract_Gaze_MO(excelLocation, dataLocation, ptID, savePreProcLocation , stimuliLOC);
clc

%% STEP 5 Run eyeTrackProc funciton ------ DONE %% Originally commented out
% STEP 5: Run eyeTRACKproc.m f(x) 

savePreProcLocation = [dataLocation , filesep ,ptID , filesep ,...
    'Eye-tracking\Processed\eyeDATA'];

saveCleanLocation = [dataLocation , filesep ,ptID , filesep ,...
    'Eye-tracking\Processed\eyeDATA\cleaned_eyeDATA'];

eyeTRACKproc_PupilSize_MO(saveCleanLocation, savePreProcLocation, ptID);

%% STEP 5.1 Run GAZE eyeTrackProc funciton ------ CURRENT STEP ON 10/4/2024

savePreProcLocation = [dataLocation , filesep ,ptID , filesep ,...
    'Eye-tracking\Processed\eyeDATA'];

saveCleanLocation = [dataLocation , filesep ,ptID , filesep ,...
    'Eye-tracking\Processed\eyeDATA\cleaned_eyeDATA'];

eyeTRACKproc_PupilLocation_MO(saveCleanLocation, savePreProcLocation, ptID);