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

lfreq = 0.5;

%% Test folder
folderDir = 'Z:\Tyner_K_Projects\MEM_ORDER\DataFolder\MW19\'; % 30, 27, 26, 25, 24, 23, 22, 21, 20, 19

%% Construct object
obj = memOrder(folderDir);

%% Select channels
selectChannels(obj,region,hemi);

%% Identify bad channels
identifyBads(obj);

%% High pass filter
dataFilter(obj,lfreq);

%% Bipolar referencing
bipolarReference(obj);

%% Identify TTLs
identifyTTLs(obj);

%% Identify the boundaries for the individual trials
identifyBoundaries(obj);

%% Break the data up into sections and epochs
allocateData(obj);

%% Baseline correction
baselineCorrection(obj);

%% Plot LFPs for each task
%plotPhases(obj);

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