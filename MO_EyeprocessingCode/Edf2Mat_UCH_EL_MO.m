% ---Function 'Edf2Mat_UCH'---
% 
% Convert .edf file to .mat file 
% Uses Edf2Mat Matlab toolbox found at github.com/mzhuang666/edf-converter
%      NOTE: Need to have this toolbox downloaded and on path 
%
% Inputs: 1) file to convert 2) Patient ID, 3) block ('en','sr','td'), 
% 4) variant ('1', '2', or '3') 
%
% Output: File saved with patient ID, block type, and date of recording
%
% Example function call: Edf2Mat_UCH('NO20221615110.edf', 'MW9', 'en' , directory1,...
%   directory2)
%
% Initial version: May 11 2022
% Update 1: July 27 2022
% Update 2: August 2 2022
% Update 3: August 26 2024
%
%
% Marielle L. Darwin & John A. Thompson 

function [eyeProcName] = Edf2Mat_UCH_EL_MO(edfFile, patientID, block, fileDIR, saveDIR)
% edfFile = 'NO20221615110.edf';
% patientID = 'MW9';
% block = 'en';


% Set path structure
basePath = fileDIR;
savePath = saveDIR;

cd(basePath);

% Construct new file name
% eyeSp = split(edfFile,'.');
% eyeRM = extractAfter(eyeSp(1),'NO');
filename = [strcat('eyetrackMO_', patientID, '_', block, '.mat')];

% Convert chosen .edf file to .mat file and save
%edfRAW = Edf2Mat(edfFile);
edfRAW = edfmex(edfFile);
cd(savePath)

% Convert edfRAW to struct
edfRawStruct = edfRAW;

save(filename,'edfRawStruct');
eyeProcName = filename;
end