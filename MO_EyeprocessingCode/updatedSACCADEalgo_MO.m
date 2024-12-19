function [] = updatedSACCADEalgo_MO(mainLOC, ptID)

% Location of eye tracking folders
cd(mainLOC)
[outFOLDS] = getfiles(mainLOC,1,nan);

% Match specific patient ID to location in outFOLDS
ptIdx = strcmp(outFOLDS, ptID);

tempCASEd = [mainLOC , filesep , outFOLDS{ptIdx}];

clGAZEdataLOC = [tempCASEd , filesep , 'Eye-tracking\Processed\eyeDATA\cleaned_eyeDATA'];

cd(clGAZEdataLOC)

matFileL1 = dir('*.mat');
matFileL2 = {matFileL1.name};
matFileL3 = matFileL2{contains(matFileL2,'clgaze')};

checkSACCADE_STROOP(matFileL3)


end




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
