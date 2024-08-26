%% Clear
clear;
clc;

%% Add path
addpath('C:\Users\Kevin_Tyner\Documents\MATLAB\matnwb-2.6.0.2');

%% Set properties
freq = [4;8]; % 4;8 8;12
boundary = ["SB"]; % "SB", "HB", "NB"
channels = ["anterior hippocampus"];
hemi = ["L"];
method = ["granger"]; % "granger", "coherence"

%% Test folder
folder = 'Z:\Tyner_K_Projects\MEM_ORDER\DataFolder\MW18\NWBProcessing\NWB_Data\';

%% Constructor
obj = test_MemOrder(folder);