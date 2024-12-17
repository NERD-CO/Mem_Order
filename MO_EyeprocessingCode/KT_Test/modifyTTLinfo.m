function [TTLinfo] = modifyTTLinfo(TTLinfo,ptID,moSessionMatFname,respMat)

    %% function [TTLinfo] = modifyTTLinfo(TTLinfo,ptID,moSessionMatFname)
    % This function modifies the trials of TTLinfo to add the corresponding
    % NLX indices for the events

    %% Add path to read files
    addpath('C:\Users\Kevin_Tyner\Documents\MATLAB\matnwb-2.6.0.2\');

    %% Grab the info file
    info_file = append('Z:\Tyner_K_Projects\MEM_ORDER\DataFolder\',ptID,'\info.xlsx');
    fileIDs = readtable(info_file);

    %% Identify the file to load
    I = strcmp(fileIDs.Behav,moSessionMatFname);
    fname = fileIDs.File{I};

    %% Find the index of the logical array
    idx = find(I,1);

    %% Count the number of times the file name has occurred
    count = sum(strcmp(fileIDs.File(1:idx), fname));

    %% Create the file name to load
    fileName = append('Z:\Tyner_K_Projects\MEM_ORDER\DataFolder\',ptID,'\NWBProcessing\NWB_Data\',fname);

    %% Load the data
    tmp = nwbRead(fileName);

    %% Load macrowire timestamps
    timestamps = tmp.processing.get('ecephys').nwbdatainterface.get...
        ('LFP').electricalseries.get('MacroWireSeries').timestamps.load;

    %% Downsample timestamps
    LFP_timestamps = downsample(timestamps,8);

    %% Extract event key
    eventTimes = tmp.acquisition.get('events').timestamps.load();
    temp_eventIDs = cellstr(tmp.acquisition.get('events').data.load());

    %% Convert eventIDs from hexadecimal
    I = contains(temp_eventIDs,'TTL');
    eventTimes = eventTimes(I);
    temp_eventIDs2 = temp_eventIDs(I);
    TTL_task = extractBetween(temp_eventIDs2,'(',')');
    eventIDs = cellfun(@(x) hex2dec(x),TTL_task,'UniformOutput',true);

    %% Find the indices
    LFPidx = zeros(length(eventTimes),1);
    for ii = 1:length(eventTimes)
        [~,LFPidx(ii,1)] = min(abs(eventTimes(ii,1) - LFP_timestamps));
    end

    %% Find corresponding start and stop indices
    x = find(eventIDs == 61); % task start
    y = find(eventIDs == 60); % task end

    %% Quick fix for start and stop if needed
    if length(y) < length(x) % No stop
        a = length(x);
        y(a,1) = length(eventIDs);
    end

    if length(x) < length(y) % No start
        x(1,1) = 1;
    end

    %% Reduce to data for individual file
    eventIDs = eventIDs(x(count,1):y(count,1),1);
    %eventTimes = eventTimes(x(count,1):y(count,1),1);
    LFPidx = LFPidx(x(count,1):y(count,1),1);

    %% Grab the TTLs
    crossOnIdx = LFPidx(eventIDs == 11);
    clipOnIdx = LFPidx(eventIDs == 1);
    clipOffIdx = LFPidx(eventIDs == 2);
    questionIdx = LFPidx(eventIDs == 3);

    %% Remove select trials
    if strcmp(ptID,'MW27') && contains(moSessionMatFname,'timeDiscrim')
        questionIdx(23) = [];
    end

    %% Get response times and convert to samples
    responseIdx = NaN(length(respMat),1);
    for ii = 1:length(respMat)
        if ~isempty(respMat(ii).respTime)
            responseIdx(ii,1) = respMat(ii).respTime - respMat(ii).QuesStart;
        else
            continue
        end
    end
    responseIdx = questionIdx + floor(responseIdx .* 500); % 500 Hz sampling rate

    %% Loop to add NLX info
    for ii = 1:length(TTLinfo)
        if istable(TTLinfo{ii,1})
            TTLinfo{ii,1}{2,4} = crossOnIdx(ii,1);
            TTLinfo{ii,1}{3,4} = clipOnIdx(ii,1);
            TTLinfo{ii,1}{4,4} = clipOffIdx(ii,1);
            TTLinfo{ii,1}{5,4} = questionIdx(ii,1);
            TTLinfo{ii,1}{6,4} = responseIdx(ii,1);
        end
    end

end