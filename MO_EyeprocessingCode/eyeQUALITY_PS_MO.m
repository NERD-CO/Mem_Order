function [] = eyeQUALITY_PS_MO(cleanedDataLOC, ptID)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

close all
% [For every file in eyeDATA, loop through to run analyses?]
% save cleaned files in new folder to only have raw data in this folder
% need to cd to new folder to save

% NO SHORTEN - because NOW Table INCLUDES time MARKERS

% CD to mainPath for raw data
cd(cleanedDataLOC);

% % Get contents of eyeDATA as a variable
% test = dir('*.mat');
% fileList = {test.name};

% for loop to loop throguh fileList
%for fileS = 1:length(fileList)
eyeData_pt = append('cl_eyeData_', ptID,'.mat');
tempFile_name = eyeData_pt;

% Load in file
load(tempFile_name, 'outInfo');

% Set tempEye- loop through every eye in every variant within variantS
% Go into variantS and determine # variants there
sessSnum = length(fieldnames(outInfo));

% Other option: Extract field names of variantS ahead of time
sessSfieldN = fieldnames(outInfo);
% Extract field names of each variant in variantS

for i = 1:sessSnum

    % name of variant
    curSession = outInfo.(sessSfieldN{i});

    for eyE = 1:4

        % FOR ERROR CHECKING
        % disp(eyE)

        switch eyE
            case 1
                titleINFO = [sessSfieldN{i} , ' Left Eye'];
                plotQLE(curSession.leftEYE.oT_pupilS_rawCL,[1 0 0],titleINFO);
            case 2
                titleINFO = [sessSfieldN{i} , ' Right Eye'];
                plotQLE(curSession.rightEYE.oT_pupilS_rawCL,[0 0 0],titleINFO);

        end
    end


end





end




function [] = plotQLE(inDATA , inCOLOR,titLEE)

figure;
for nfi = 1:length(inDATA)
    tmpEFi = inDATA{nfi};
    plot(tmpEFi,'Color',[inCOLOR 0.3],'LineWidth',1.5)
    hold on
end
title(titLEE)

end