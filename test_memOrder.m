%% Clear
clear;
clc;

%% Add path
addpath('C:\Users\Kevin_Tyner\Documents\MATLAB\matnwb-2.6.0.2');

%% Set properties
freq = [4;8]; % 4;8 8;12
boundary = ["SB"]; % "SB", "HB", "NB"
region = ["anterior hippocampus"];
hemi = ["L"];
method = ["granger"]; % "granger", "coherence"

lfreq = 1; % 0.5
%hfreq = 250;

%% Test folder
folderDir = 'Z:\Tyner_K_Projects\MEM_ORDER\DataFolder\MW30\'; % 30, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18

%% Construct object
obj = memOrder(folderDir);

%% Load eye data
loadEyeData(obj);

%% Fix channel names
fixChNames(obj);

%% Identify bad channels
identifyBads(obj);

%% Bipolar reference
bipolarReference(obj);

%% Filter
dataFilter(obj,lfreq);

%% Identify TTLs
identifyTTLs(obj);

%% Identify the boundaries for the individual trials
identifyBoundaries(obj);

%% Break the data up into sections and epochs
allocateData(obj);

%% Baseline correction
baselineCorrection(obj);

%% Z-score the epochs
calculateZscore(obj);

%% Notes
% All subsequent functions are optional.

%% Plot LFPs for each task
plotPhases(obj);

%% Average spectrogram
averageSpectrogram(obj);

%% Plot spectrograms
% Only run this function if the averageSpectrogram function was run
plotSpectrograms(obj);

%% Plot frequency bands
plotFrequencyBands(obj);

%% Examine theta phase relation to saccades
thetaSaccade(obj,2,10);

%% Compute power
computePower(obj,freq);

%% Compute statistics
computeStats(obj);

%% Functional connectivity
%functionalConnectivity(obj,method);

%% Analyze significant changes in FC
%analyzeFC(obj,boundary);

%% SPOOOF
SPOOOF(obj,0.5,250);