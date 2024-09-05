function[behavFILE] = MO_TxttoMat(textloc,patientID,block, tempCASEd)
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


% Load text block by reading in every line of .txt file
inputFile=fopen(textloc);       % Open file in Matlab
tLine=fgetl(inputFile);         % Read in 1st line of .txt file
%tLine=(fgetl(inputFile)+1);    % Read in 2nd line of .txt file (where ttls start)
allLines={};                    % Create empty cell array
counter=1;                      % Create a counter
while ischar(tLine)
    allLines{counter}=tLine;     % Fill cell array with txt lines
    counter=counter+1;           % Update counter
    tLine=fgetl(inputFile);      % Move through lines in txt file
    %disp(num2str(counter))       % Display counter
end

% Create table from .txt file info
%
% Create cellfun to split elements in 'allLines' & output to a cell array
taskElements_temp=cellfun(@(x) strsplit(x,';'),allLines,'UniformOutput',false);
% Cut 1st column of taskElements
taskElements = taskElements_temp(2:width(taskElements_temp));
% Create cellfun to "loop" through 'taskElements' to extract 2nd element
TTL_ID=cellfun(@(x) x{2},taskElements,'UniformOutput',false);
% Create cellfun to "loop" through 'taskElements' to extract 1st element
timestamp=cellfun(@(x) x{1},taskElements,'UniformOutput',false);
% Interpreter to translate TTL identifiers
[TTLdescription]=TTLinterpret(TTL_ID);

% Create output table
taskData=table(transpose(TTL_ID),transpose(TTLdescription),transpose(timestamp),...
    'VariableNames',{'TTLvalue','TTL_ID','Timestamp'});

% Store in .mat file
outData.taskinformation = taskData;
outData.patientID = patientID;
outData.moSession = block;

% Save file
cd(tempCASEd);
savefilename=[block,'_',patientID,'.mat'];
save(savefilename,'outData');
behavFILE = savefilename;
end

% Create subfunction to make an interpreter for TTL values
function[TTLdescription]=TTLinterpret(TTLnumber)
%Create table for TTL value key
TTL_ID={'TaskOnset','FixCross','ClipOnset','ClipOffset','Probe','Response','TaskOffset','AfterTTLdelay'};
TTLvalue={'61','11','1','2','3','4','60','0'};

keyTable=table(transpose(TTL_ID),transpose(TTLvalue),'VariableNames',{'TTLI','TTLN'});

% Create empty array to preallocate output
TTLdescription=cell(size(TTLnumber));

% Loop through each TTL number in input argument
for i=1:length(TTLnumber)
    tempN=TTLnumber{i};
    % Find matching ID for TTL value from keyTable

    % KT add - 8/27/24
    if ~ismember(tempN,keyTable.TTLN)
        TTLdescription{i} = 'NaN';
    else
        TTLlabel=keyTable.TTLI(ismember(keyTable.TTLN,tempN));
        TTLdescription{i}=TTLlabel{1};
    end

    %TTLlabel=keyTable.TTLI(ismember(keyTable.TTLN,tempN));
    %TTLdescription{i}=TTLlabel{1}; %1=value is true, 0=value is false
end

end