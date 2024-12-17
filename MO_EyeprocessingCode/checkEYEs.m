%% Clear
clear
clc

%% Subject Location
pID = 'MW30';
subLoc = append('Z:\Tyner_K_Projects\MEM_ORDER\DataFolder\',pID,'\Eye-tracking\Processed\eyeDATA\cleaned_eyeDATA');

%% Check eye quality
eyeQUALITY_PS_MO(subLoc, pID)