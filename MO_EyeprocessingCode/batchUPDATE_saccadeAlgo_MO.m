cd('E:\Dropbox\NEWOLD_Delay_2025\SUBJECT_Data\')
subDIRall1 = dir();
subDIRall2 = {subDIRall1.name};
subDIRall3 = subDIRall2(~ismember(subDIRall2,{'.','..'}));


for ssii = 1:length(subDIRall3)

    subTEMPi = subDIRall3{ssii};

    addGAZEv2(subTEMPi)

end

%%
mainLOC = 'E:\Dropbox\StroopEMUTest_2024';
ptID = 'MW40';
updatedSACCADEalgo_STROOP(mainLOC, ptID)