function[behavFILE] = MO_TxttoMat(textloc,patientID,block, tempCASEd)

%% function[behavFILE] = MO_TxttoMat(textloc,patientID,block, tempCASEd)
% Convert .txt file from behavioral session into .mat file
% Outputs table to view TTL values and timestamps
%
% John A. Thompson & Marielle L. Darwin
% Created July 7 2021 | Updated April 24 2023
%
%
% Events description within .txt file
% The events coorespond to the TTL markers for each trial
% For the learning trials, the TTL markers are the following:
% 55 = start of the experiment, 1 = stimulus ON, 2 = stimulus OFF, 3 = Question Screen Onset,
% 20 = Yes, 21 = NO, 6 = End of Delay after Response, 66 = End of Experiment.
% For the recognition trials, the TTL markers are the following:
% 55 = start of the experiment, 1 = stimulus ON, 2 = stimulus OFF, 3 = Question Screen Onset,
% 31:36 = Confidence in response, 66 = End of Experiment;
%
% Key for TTL_ID shorthand in the 'TTLinterpret' subfunction:
% 'NS'=no, sure; 'NLS'=no, less sure; 'NVU'=no, very unsure
% 'YS'=yes, sure; 'YLS'=yes, less sure; 'YVU'=yes, very unsure

%% Load text block by reading in every line of .txt file
inputFile=fopen(textloc);
tLine=fgetl(inputFile);
allLines={};
counter = 1;
while ischar(tLine)
    allLines{counter}=tLine;
    counter=counter+1; 
    tLine=fgetl(inputFile);
end

%% Create table from .txt file info
taskElements_temp=cellfun(@(x) strsplit(x,';'),allLines,'UniformOutput',false);
taskElements = taskElements_temp(2:width(taskElements_temp));
TTL_ID=cellfun(@(x) x{2},taskElements,'UniformOutput',false);
timestamp=cellfun(@(x) x{1},taskElements,'UniformOutput',false);
[TTLdescription]=TTLinterpret(TTL_ID);

%% Create output table
taskData=table(transpose(TTL_ID),transpose(TTLdescription),transpose(timestamp),...
    'VariableNames',{'TTLvalue','TTL_ID','Timestamp'});

%% Store in .mat file
outData.taskinformation = taskData;
outData.patientID = patientID;
outData.moSession = block;

%% Save file
cd(tempCASEd);
savefilename=[block,'_',patientID,'.mat'];
save(savefilename,'outData');
behavFILE = savefilename;

end