function [] = Extract_Eye_EDF_MO(excelLOC , mainLOC, ptID)
% The events coorespond to the TTL markers for each trial::::::
%  61 = Task onset
%  60 = Task offset
%  11 = Fixation cross
%  1  = Clip Onset
%  2  = Clip Offset
%  3  = Probe
%  4  = Response
%  0  = afterTTLdelay
%  Differences between 
%  MO Session: Encoding
%            : Scene Recognition
%            : Time Discrimination
% Updated : August 26 2024

% Location of variant.xlsx
cd(excelLOC)
varTable = readtable('Eye_file_information.xlsx');
% Location of eye tracking folders
cd(mainLOC)
[outFOLDS] = getfiles(mainLOC,1,nan);

% Match specific patient ID to location in outFOLDS
ptIdx = find(strcmp(outFOLDS, ptID));

rawEYEdir = [mainLOC , filesep , outFOLDS{ptIdx}, filesep ,...
    'Eye-tracking\Raw'];

saveEYEdir = [mainLOC , filesep , outFOLDS{ptIdx} , filesep ,...
    'Eye-tracking\Processed'];
idTab = varTable(matches(varTable.Subject,outFOLDS{ptIdx}),:);
cd(rawEYEdir)

avail_MOsessions = ~matches(idTab.EyeData,'NaN');
all_MOsessions = idTab.TaskID(avail_MOsessions);
all_MOfiles = idTab.EyeData(avail_MOsessions);

for mo_i = 1:length(all_MOsessions)

    % CREATE the RAW EDF file from eye data
    Edf2Mat_UCH_EL_MO(all_MOfiles{mo_i}, outFOLDS{ptIdx}, all_MOsessions{mo_i},...
        rawEYEdir, saveEYEdir);

end


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


