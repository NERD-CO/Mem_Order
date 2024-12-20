cd('E:\Dropbox\NEWOLD_Delay_2025\SUBJECT_Data\')
subDIRall1 = dir();
subDIRall2 = {subDIRall1.name};
subDIRall3 = subDIRall2(~ismember(subDIRall2,{'.','..'}));


for ssii = 1:length(subDIRall3)

    subTEMPi = subDIRall3{ssii};

    addGAZEv2(subTEMPi)

end

%%
mainLOC = 'Z:\Tyner_K_Projects\MEM_ORDER\DataFolder';
ptID = 'MW18';
updatedSACCADEalgo_MO(mainLOC, ptID)


%%
cd(mainLOC)
foldDIR = dir();
foldDIR2 = {foldDIR.name};
subLIST = foldDIR2(~ismember(foldDIR2,{'.','..'}));
subLIST = subLIST(2:end);

for subII = 1:length(subLIST)

    ptID = subLIST{subII};
    updatedSACCADEalgo_MO(mainLOC, ptID)
    disp([ptID , ' done'])

end

