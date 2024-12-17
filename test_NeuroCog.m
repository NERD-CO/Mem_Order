%% Clear
clear;
clc;

%% Add path
addpath('C:\Users\Kevin_Tyner\Documents\MATLAB\matnwb-2.6.0.2');

%% Properties
refScheme = 'bipolar'; % can be 'bipolar' or 'laplace'
lfreq = 1;
hfreq = 200;
type = 'fir'; % can be 'fir' or 'iir'
order = 256; % for 'fir', order = 256; for 'iir', order = 2

%% Test folder
folderDir = 'Z:\Tyner_K_Projects\NEW_OLD_DELAY\SUBJECT_Data\MW12_PY22NO12\Variant3\';
%folderDir = 'Z:\Tyner_K_Projects\MEM_ORDER\DataFolder\MW23\'; % 19, 21, 22, 23

%% Constructor object
obj = NeuroCog(folderDir);

%% Load eye data
loadEyeData(obj);

%% Identify bad channels
identifyBads(obj);

%% Reference
dataReference(obj,refScheme);

%% Filter
dataFilter(obj,lfreq,hfreq,type,order);

%% Identify TTLs
identifyTTLs(obj);

%% Identify boundaries (MemOrder only)
identifyBoundaries(obj);

%% Allocate the data
allocateData(obj);

%% Baseline correction and Z-scoring
baselineCorrectZscore(obj);

%% Plot the task phases - pure sanity check
%plotPhases(obj);

%% Theta saccades
saccadeAnalysis(obj,2,10,type,order);