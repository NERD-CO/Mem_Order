
% STEP 1: Change path directories. Choose the NAS case but update the paths
% the first time you use it to make sure they are correct.

getPCname = getenv('COMPUTERNAME');

switch getPCname
    case 'DESKTOP-FAGRV5G'
        % % Code location
        % codeLocation = 'C:\Users\Admin\Documents\Github\Neurocuanalysis\Mem_Order_analysis';
        % addpath(codeLocation)
        % % Mex file for converting eye tracking data
        % edfMexELloc = 'C:\Users\Admin\Documents\Github\Neurocuanalysis\edfmex';
        % edfCheck2 = which('edfmex.mexa64');
        % if isempty(edfCheck2)
        %     addpath(genpath(edfMexELloc));
        % end
        % 
        % % OPTIONAL - color maps
        % cbrewLOC = 'C:\Users\Admin\Documents\Github\Neurocuanalysis\cbrewerALL';
        % addpath(genpath(cbrewLOC))
        % 
        % % File with subject information 
        % excelLocation = 'Z:\Tyner_K_Projects\MEM_ORDER';
        % % Data folder location
        % dataLocation = 'Z:\Tyner_K_Projects\MEM_ORDER\DataFolder';
        % % Save locations for Eye data
        % savePreProcLocation = [dataLocation , filesep , 'eyeDATA'];
        % saveCleanLocation = [savePreProcLocation , filesep , 'cleaned_eyeDATA'];
        % cd(codeLocation)
    case 'DESKTOP-I5CPDO7'
        % Code location
        codeLocation = 'C:\Users\Admin\Documents\Github\Neurocuanalysis\Mem_Order_analysis';
        addpath(codeLocation)
        % Mex file for converting eye tracking data
        edfMexELloc = 'C:\Users\Admin\Documents\Github\Neurocuanalysis\edfmex';
        edfCheck2 = which('edfmex.mexa64');
        if isempty(edfCheck2)
            addpath(genpath(edfMexELloc));
        end

        % OPTIONAL - color maps
        cbrewLOC = 'C:\Users\Admin\Documents\Github\Neurocuanalysis\cbrewerALL';
        addpath(genpath(cbrewLOC))

        % File with subject information 
        excelLocation = 'Z:\Tyner_K_Projects\MEM_ORDER';
        % Data folder location
        dataLocation = 'Z:\Tyner_K_Projects\MEM_ORDER\DataFolder';
        % Save locations for Eye data
        savePreProcLocation = [dataLocation , filesep , 'eyeDATA'];
        saveCleanLocation = [savePreProcLocation , filesep , 'cleaned_eyeDATA'];
        cd(codeLocation)
    case 'MLD'

    case 'NAS'
 

end

%% STEP 2: Change ptID to be specific to pt 
ptID = 'MW18';

%% STEP 3 CONVERT From EDF to MAT
Extract_Eye_EDF_MO(excelLocation , dataLocation, ptID)
clc

%% STEP 4 Run Initial_EyeAnalysis function

savePreProcLocation = [dataLocation , filesep ,ptID , filesep ,...
    'Eye-tracking\Processed\eyeDATA'];

EyeAnalysis_DataExtract_MO(excelLocation, dataLocation, ptID, savePreProcLocation);
clc

%% STEP 4.1 RUN initial eye position function

EyeAnalysis_DataExtract_EL_Gaze(excelLocation, dataLocation, ptID, savePreProcLocation);
clc
