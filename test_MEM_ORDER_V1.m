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
folder = 'Z:\Tyner_K_Projects\MEM_ORDER\DataFolder\MW19\NWBProcessing\NWB_Data\'; %26
fnames = dir(folder);
fnames = natsortfiles(fnames); % MATLAB add-on
fnames = {fnames(~[fnames.isdir]).name}';
I = contains(fnames,'filter');
fnames = fnames(I);

%% Create taskID
taskID = cell(length(fnames),1);

%% Loop through to identify the task per file
for ii = 1:length(fnames)

    %% Create file name
    fname = append(folder,fnames{ii,1});

    %% Load the data
    data = nwbRead(fname);

    %% Get the TTLs
    temp_eventIDs = cellstr(data.acquisition.get('events').data.load());

    %% Convert
    J = contains(temp_eventIDs,'TTL');
    temp_eventIDs2 = temp_eventIDs(J);
    TTL_task = extractBetween(temp_eventIDs2,'(',')');
    eventIDs = cellfun(@(x) hex2dec(x),TTL_task,'UniformOutput',true);

    %% Give file a taskID
    if sum(eventIDs == 11) == 90
        if ~any(strcmp(taskID,'encoding'))
            taskID{ii,1} = 'encoding';
        else
            taskID{ii,1} = 'timeDiscrim';
        end
    elseif sum(eventIDs == 11) == 180
        taskID{ii,1} = 'sceneRecog';
    end

end

%% Loop through and preprocess individual files
for ii = 1:length(fnames)

    %% Check if there is a file tag
    if ~isempty(taskID{ii,1})

        %% Create the file name
        fname = append(folder,fnames{ii,1});

        %% Construct the object
        obj = MEM_ORDER(fname,taskID{ii,1});

        %% Select channels
        selectChannels(obj,channels,hemi);

        %% Identify bad channels
        identifyBads(obj);

        %% Perform bipolar referencing
        bipolarMontage(obj);

        %% Identify where TTLs occur in LFP
        identifyTTLs(obj);

        %% Identify task boundaries
        identifyBoundaries(obj);

        %% Allocate the data
        allocateData(obj);

        %% Baseline correction
        baselineCorrection(obj);

        %% Plot the LFP
        plotPhases(obj);

        %% Compute power
        computePower(obj,freq);

        %% Statistical analysis
        computeStats(obj,boundary);

        %% Functional Connectivity
        functionalConnectivity(obj,method);

        %% Analyze functional connectivity
        analyzeFC(obj);
        x = 0;


    end
end