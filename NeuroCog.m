classdef NeuroCog < handle

    %% Properties
    properties(SetAccess = private)
        subID
        folderID
        task
        data
    end

    %% Methods
    methods

        %% Constructor
        function obj = NeuroCog(folderDir)

            %% Notes
            % This constructor object differentially loads in data for the
            % New-Old delay task and the memOrder task.

            %% Set folderID and subID
            if contains(folderDir,'SUBJECT_Data')
                obj.subID = extractBefore(extractAfter(folderDir,'SUBJECT_Data\'),'\V');
                obj.task = extractAfter(extractBefore(folderDir,'\SUBJECT_Data'),'Projects\');
            else
                obj.subID = extractBefore(extractAfter(folderDir,'DataFolder\'),'\');
                obj.task = extractAfter(extractBefore(folderDir,'\DataFolder'),'Tyner_K_Projects\');
            end
            obj.folderID = folderDir;

            %% Load Pt summary file for Mem Order
            if ~contains(folderDir,'SUBJECT_Data')
                pt_info = dir(fullfile(folderDir,'*.xlsx'));
                pt_info = append(folderDir,pt_info.name);
                fileIDs = readtable(pt_info);

                %% Identify the tasks the patient performed
                info = fileIDs.Properties.VariableNames;
                count = 1;
                tasks = cell(3,1);
                if strcmp(obj.task,'MEM_ORDER')
                    for ii = 2:4
                        if any(fileIDs{:,ii} == 1)
                            tasks{count,1} = info{1,ii};
                            count = count + 1;
                        end
                    end
                end
                tasks = tasks(~cellfun('isempty',tasks));

                %% Get file names
                fileNames = fileIDs.File;
                behavNames = fileIDs.Behav;
                eyeNames = fileIDs.Eye_file;
                eye = fileIDs.Eye;

            end

            %% Load the data
            if strcmp(obj.task,'NEW_OLD_DELAY')

                %% Find the mat files
                files = dir(fullfile(obj.folderID,'*.mat'));
                names = {files.name}';
                obj.data.Learn.meta_file = append(obj.folderID,string(names{contains(names,'Learn')}));
                obj.data.Recog.meta_file = append(obj.folderID,string(names{contains(names,'Recog')}));

                %% Set the tasks
                tasks = {"Learn";"Recog"};

                %% Loop over the tasks
                for ii = 1:length(tasks)

                    %% Load the metadata
                    obj.data.(tasks{ii,1}).info = load(obj.data.(tasks{ii,1}).meta_file,"sessionINFO");
                    obj.data.(tasks{ii,1}).info = obj.data.(tasks{ii,1}).info.sessionINFO;

                    %% Create NWB file
                    nwbName = append(obj.folderID,obj.data.(tasks{ii,1}).info.NWBfile);

                    %% Print
                    fprintf('Ephys data loaded for sub %s and task %s..\n',obj.subID,tasks{ii,1});

                    %% Load NWB file
                    tmp = nwbRead(nwbName);

                    %% Load macrowire timestamps
                    timestamps = tmp.processing.get('ecephys').nwbdatainterface.get...
                        ('LFP').electricalseries.get('MacroWireSeries').timestamps.load;

                    %% Downsample timestamps
                    obj.data.(tasks{ii,1}).LFP_info.timestamps = downsample(timestamps,8);

                    %% Convert time stamps
                    obj.data.(tasks{ii,1}).LFP_info.converted_timestamps = (obj.data.(tasks{ii,1}).LFP_info.timestamps - ...
                        obj.data.(tasks{ii,1}).LFP_info.timestamps(1,1))/(1e6);

                    %% Get the sampling frequency
                    LFP_sessionInfo = tmp.processing.get('ecephys').nwbdatainterface.get...
                        ('LFP').electricalseries.get('MacroWireSeries');
                    obj.data.(tasks{ii,1}).fs = str2double(cell2mat(extractBetween(LFP_sessionInfo.description,'= ',':')));

                    %% Get voltages for macrowires and convert
                    obj.data.(tasks{ii,1}).LFP_info.data = (double(tmp.processing.get('ecephys').nwbdatainterface.get...
                        ('LFP').electricalseries.get('MacroWireSeries').data.load)) .* LFP_sessionInfo.data_conversion;

                    %% Channel IDs
                    chanLabels = cellstr(tmp.general_extracellular_ephys_electrodes.vectordata.get('label').data.load()); %use MA only!
                    MAchan = find(contains(chanLabels,'MA_'));
                    chanID = cellstr(tmp.general_extracellular_ephys_electrodes.vectordata.get('location').data.load());
                    hemisphere = cellstr(tmp.general_extracellular_ephys_electrodes.vectordata.get('hemisph').data.load());
                    shortBnames = cellstr(tmp.general_extracellular_ephys_electrodes.vectordata.get('shortBAn').data.load());
                    wireID = tmp.general_extracellular_ephys_electrodes.vectordata.get('channID').data.load();

                    obj.data.(tasks{ii,1}).LFP_info.chanID = chanID(MAchan);
                    obj.data.(tasks{ii,1}).LFP_info.chanHemi = hemisphere(MAchan);
                    obj.data.(tasks{ii,1}).LFP_info.chanSname = shortBnames(MAchan);
                    obj.data.(tasks{ii,1}).LFP_info.wireID = wireID(MAchan);

                    %% Find unique brain regions and turn into a table
                    tmpBRegUni = unique(obj.data.(tasks{ii,1}).LFP_info.chanSname);
                    hemiTemp = cell(length(tmpBRegUni),1);
                    longName = cell(length(tmpBRegUni),1);
                    for ui = 1:length(tmpBRegUni)
                        tmpIND = find(matches(obj.data.(tasks{ii,1}).LFP_info.chanSname,tmpBRegUni{ui}),1,'first');
                        hemiTemp{ui} = obj.data.(tasks{ii,1}).LFP_info.chanHemi{tmpIND};
                        longName{ui} = obj.data.(tasks{ii,1}).LFP_info.chanID{tmpIND};
                    end
                    obj.data.(tasks{ii,1}).LFP_info.brTABLE = table(tmpBRegUni,hemiTemp,longName,'VariableNames',{'SEEGele',...
                        'Hemisphere','LongBRname'});

                    %% Assign LFP to each brain region
                    obj.data.(tasks{ii,1}).LFP_info.region_LFP = cell(length(obj.data.(tasks{ii,1}).LFP_info.brTABLE.SEEGele),1);
                    for k = 1:height(obj.data.(tasks{ii,1}).LFP_info.brTABLE)
                        reg = strcmp(obj.data.(tasks{ii,1}).LFP_info.brTABLE.SEEGele{k},obj.data.(tasks{ii,1}).LFP_info.chanSname);
                        obj.data.(tasks{ii,1}).LFP_info.region_LFP{k,1} = obj.data.(tasks{ii,1}).LFP_info.data(reg,:);
                        obj.data.(tasks{ii,1}).LFP_info.regionChan{k,1} = obj.data.(tasks{ii,1}).LFP_info.chanSname(reg,1);
                    end

                    %% Fix channel names
                    for jj = 1:length(obj.data.(tasks{ii,1}).LFP_info.regionChan)
                        for kk = 1:length(obj.data.(tasks{ii,1}).LFP_info.regionChan{jj,1})
                            obj.data.(tasks{ii,1}).LFP_info.regionChan{jj,1}(kk,1) = append(obj.data.(tasks{ii,1}).LFP_info.regionChan{jj,1}(kk,1),num2str(kk));
                        end
                    end

                end

            elseif strcmp(obj.task,'MEM_ORDER')

                %% Preallocate
                count = 1;
                matchedNames = cell(length(tasks),1);

                %% Loop over files
                for ii = 1:length(fileNames)

                    %% Check to see if file has been used before
                    if any(strcmp(fileNames{ii,1},matchedNames))
                        count = count + 1;
                    else
                        count = 1;
                    end

                    %% Add file name to check
                    matchedNames{ii,1} = fileNames{ii,1};

                    %% Add file names to obj
                    obj.data.(tasks{ii,1}).nwbName = append(folderDir,'NWBProcessing\NWB_Data\',fileNames{ii,1});
                    obj.data.(tasks{ii,1}).behavName = append(folderDir,'Behavioral_Data\Raw\',behavNames{ii,1});
                    obj.data.(tasks{ii,1}).eyeName = append(folderDir,'Eye-tracking\Processed\eyeDATA\cleaned_eyeDATA\',eyeNames{ii,1});
                    obj.data.(tasks{ii,1}).eye2use = eye{ii,1};

                    %% Print
                    fprintf('Ephys data loaded for sub %s and task %s..\n',obj.subID,tasks{ii,1});

                    %% Load the NWB file
                    tmp = nwbRead(obj.data.(tasks{ii,1}).nwbName);

                    %% Load behavioral data
                    tmpStruct = load(obj.data.(tasks{ii,1}).behavName);
                    tmpNames = fieldnames(tmpStruct);
                    J = contains(tmpNames,'respMat');
                    tmpStruct.newrespMat = tmpStruct.(tmpNames{J,1});
                    tmpStruct = rmfield(tmpStruct,tmpNames{J,1});
                    tmpStruct.respMat = tmpStruct.newrespMat;
                    tmpStruct = rmfield(tmpStruct,'newrespMat');
                    obj.data.(tasks{ii,1}).behav = tmpStruct; %%%%%% Need to fix respMat values
                    clear tmpStruct;

                    %% Load macrowire timestamps
                    timestamps = tmp.processing.get('ecephys').nwbdatainterface.get...
                        ('LFP').electricalseries.get('MacroWireSeries').timestamps.load;

                    %% Downsample timestamps
                    obj.data.(tasks{ii,1}).LFP_info.timestamps = downsample(timestamps,8);

                    %% Convert time stamps
                    obj.data.(tasks{ii,1}).LFP_info.converted_timestamps = (obj.data.(tasks{ii,1}).LFP_info.timestamps - ...
                        obj.data.(tasks{ii,1}).LFP_info.timestamps(1,1))/(1e6);

                    %% Get the sampling frequency
                    LFP_sessionInfo = tmp.processing.get('ecephys').nwbdatainterface.get...
                        ('LFP').electricalseries.get('MacroWireSeries');
                    obj.data.(tasks{ii,1}).fs = str2double(cell2mat(extractBetween(LFP_sessionInfo.description,'= ',':')));

                    %% Get voltages for macrowires and convert
                    obj.data.(tasks{ii,1}).LFP_info.data = (double(tmp.processing.get('ecephys').nwbdatainterface.get...
                        ('LFP').electricalseries.get('MacroWireSeries').data.load)) .* LFP_sessionInfo.data_conversion;

                    %% Channel IDs
                    chanLabels = cellstr(tmp.general_extracellular_ephys_electrodes.vectordata.get('label').data.load()); %use MA only!
                    MAchan = find(contains(chanLabels,'MA_'));
                    chanID = cellstr(tmp.general_extracellular_ephys_electrodes.vectordata.get('location').data.load());
                    hemisphere = cellstr(tmp.general_extracellular_ephys_electrodes.vectordata.get('hemisph').data.load());
                    shortBnames = cellstr(tmp.general_extracellular_ephys_electrodes.vectordata.get('shortBAn').data.load());
                    wireID = tmp.general_extracellular_ephys_electrodes.vectordata.get('channID').data.load();

                    obj.data.(tasks{ii,1}).LFP_info.chanID = chanID(MAchan);
                    obj.data.(tasks{ii,1}).LFP_info.chanHemi = hemisphere(MAchan);
                    obj.data.(tasks{ii,1}).LFP_info.chanSname = shortBnames(MAchan);
                    obj.data.(tasks{ii,1}).LFP_info.wireID = wireID(MAchan);

                    %% Find unique brain regions and turn into a table
                    tmpBRegUni = unique(obj.data.(tasks{ii,1}).LFP_info.chanSname);
                    hemiTemp = cell(length(tmpBRegUni),1);
                    longName = cell(length(tmpBRegUni),1);
                    for ui = 1:length(tmpBRegUni)
                        tmpIND = find(matches(obj.data.(tasks{ii,1}).LFP_info.chanSname,tmpBRegUni{ui}),1,'first');
                        hemiTemp{ui} = obj.data.(tasks{ii,1}).LFP_info.chanHemi{tmpIND};
                        longName{ui} = obj.data.(tasks{ii,1}).LFP_info.chanID{tmpIND};
                    end
                    obj.data.(tasks{ii,1}).LFP_info.brTABLE = table(tmpBRegUni,hemiTemp,longName,'VariableNames',{'SEEGele',...
                        'Hemisphere','LongBRname'});

                    %% Modify channel names
                    for k = 1:height(obj.data.(tasks{ii,1}).LFP_info.brTABLE)
                        reg = strcmp(obj.data.(tasks{ii,1}).LFP_info.brTABLE.SEEGele{k},obj.data.(tasks{ii,1}).LFP_info.chanSname);
                        obj.data.(tasks{ii,1}).LFP_info.regionChan{k,1} = obj.data.(tasks{ii,1}).LFP_info.chanSname(reg,1);
                    end

                    %% Fix channel names
                    for jj = 1:length(obj.data.(tasks{ii,1}).LFP_info.regionChan)
                        for kk = 1:length(obj.data.(tasks{ii,1}).LFP_info.regionChan{jj,1})
                            obj.data.(tasks{ii,1}).LFP_info.regionChan{jj,1}(kk,1) = append(obj.data.(tasks{ii,1}).LFP_info.regionChan{jj,1}(kk,1),num2str(kk));
                        end
                    end

                    %% Extract event key
                    obj.data.(tasks{ii,1}).LFP_info.eventTimes = tmp.acquisition.get('events').timestamps.load();
                    temp_eventIDs = cellstr(tmp.acquisition.get('events').data.load());

                    %% Convert eventIDs from hexadecimal
                    I = contains(temp_eventIDs,'TTL');
                    obj.data.(tasks{ii,1}).LFP_info.eventTimes = obj.data.(tasks{ii,1}).LFP_info.eventTimes(I);
                    temp_eventIDs2 = temp_eventIDs(I);
                    TTL_task = extractBetween(temp_eventIDs2,'(',')');
                    obj.data.(tasks{ii,1}).LFP_info.eventIDs = cellfun(@(x) hex2dec(x),TTL_task,'UniformOutput',true);

                    %% Find task start and stop indices
                    x = find(obj.data.(tasks{ii,1}).LFP_info.eventIDs == 61); % task start
                    y = find(obj.data.(tasks{ii,1}).LFP_info.eventIDs == 60); % task end

                    %% Find the time in the recording that corresponds to the indices
                    for jj = 1:length(obj.data.(tasks{ii,1}).LFP_info.eventTimes)
                        [~,obj.data.(tasks{ii,1}).LFP_info.LFPIdx(jj,1)] = min(abs(obj.data.(tasks{ii,1}).LFP_info.eventTimes(jj,1) - obj.data.(tasks{ii,1}).LFP_info.timestamps));
                    end

                    %% Little fix for uneven indices if needed
                    if length(y) < length(x) % if no stop
                        a = length(x);
                        y(a,1) = length(obj.data.(tasks{ii,1}).LFP_info.eventIDs);
                    end

                    if length(x) < length(y) % if no start
                        x(1,1) = 1;
                    end

                    %% Reduce indices to those relevant to the task
                    obj.data.(tasks{ii,1}).LFP_info.eventIDs = obj.data.(tasks{ii,1}).LFP_info.eventIDs(x(count,1):y(count,1),1);
                    obj.data.(tasks{ii,1}).LFP_info.eventTimes = obj.data.(tasks{ii,1}).LFP_info.eventTimes(x(count,1):y(count,1),1);
                    obj.data.(tasks{ii,1}).LFP_info.LFPIdx = obj.data.(tasks{ii,1}).LFP_info.LFPIdx(x(count,1):y(count,1),1);

                    %% Reduce data to that of the individual task
                    obj.data.(tasks{ii,1}).LFP_info.timestamps = obj.data.(tasks{ii,1}).LFP_info.timestamps(obj.data.(tasks{ii,1}).LFP_info.LFPIdx(1,1):obj.data.(tasks{ii,1}).LFP_info.LFPIdx(end,1),:);
                    obj.data.(tasks{ii,1}).LFP_info.converted_timestamps = obj.data.(tasks{ii,1}).LFP_info.converted_timestamps(obj.data.(tasks{ii,1}).LFP_info.LFPIdx(1,1):obj.data.(tasks{ii,1}).LFP_info.LFPIdx(end,1),:);
                    obj.data.(tasks{ii,1}).LFP_info.data = obj.data.(tasks{ii,1}).LFP_info.data(:,obj.data.(tasks{ii,1}).LFP_info.LFPIdx(1,1):obj.data.(tasks{ii,1}).LFP_info.LFPIdx(end,1));

                    %% Assign LFP to each brain region
                    obj.data.(tasks{ii,1}).LFP_info.region_LFP = cell(length(obj.data.(tasks{ii,1}).LFP_info.brTABLE.SEEGele),1);
                    for k = 1:height(obj.data.(tasks{ii,1}).LFP_info.brTABLE)
                        reg = strcmp(obj.data.(tasks{ii,1}).LFP_info.brTABLE.SEEGele{k},obj.data.(tasks{ii,1}).LFP_info.chanSname);
                        obj.data.(tasks{ii,1}).LFP_info.region_LFP{k,1} = obj.data.(tasks{ii,1}).LFP_info.data(reg,:);
                    end

                    %% Calculate correction factor
                    obj.data.(tasks{ii,1}).LFP_info.correction_factor = obj.data.(tasks{ii,1}).LFP_info.LFPIdx(1,1) - 1;
                    obj.data.(tasks{ii,1}).LFP_info.LFPIdx = obj.data.(tasks{ii,1}).LFP_info.LFPIdx - obj.data.(tasks{ii,1}).LFP_info.correction_factor;

                end
            end
        end

        %% loadEyeData
        function loadEyeData(obj)

            %% function loadEyeData(obj)
            % This function loads in the eye data specifically for the
            % memOrder task.

            %% Print
            fprintf('Running function loadEyeData..\n')

            %% Check if memOrder
            if strcmp(obj.task,'MEM_ORDER')

                %% If yes, create the task names
                tasks = fieldnames(obj.data);

                %% Loop over tasks
                for ii = 1:length(tasks)

                    %% Load the eye data
                    outInfo = load(obj.data.(tasks{ii,1}).eyeName);
                    outInfo = outInfo.outInfo;

                    %% Allocate to task
                    if strcmp(tasks{ii,1},'Encode')
                        obj.data.(tasks{ii,1}).eyeInfo = outInfo.encoding;
                    elseif strcmp(tasks{ii,1},'SceneRecognition')
                        obj.data.(tasks{ii,1}).eyeInfo = outInfo.sceneRecog;
                    elseif strcmp(tasks{ii,1},'TimeDiscrimination')
                        obj.data.(tasks{ii,1}).eyeInfo = outInfo.timeDiscrim;
                    end

                end
            end
        end

        %% Identify bad channels
        function identifyBads(obj)

            %% function identifyBads(obj)
            % This function identifies bad channels based on whether or not
            % their amplitudes cross a certain threshold for a specified
            % number of samples.

            %% Print
            fprintf('Running function identifyBads..\n');

            %% Get the task names
            tasks = fieldnames(obj.data);

            %% Loop over tasks to identify bads
            for ii = 1:length(tasks)

                %% Loop over wires
                for jj = 1:length(obj.data.(tasks{ii,1}).LFP_info.region_LFP)

                    %% Loop over channels
                    for kk = 1:size(obj.data.(tasks{ii,1}).LFP_info.region_LFP{jj,1},1)

                        %% Calculate channel mean and standard deviation
                        ch_mean = mean(obj.data.(tasks{ii,1}).LFP_info.region_LFP{jj,1}(kk,:));
                        ch_std = std(obj.data.(tasks{ii,1}).LFP_info.region_LFP{jj,1}(kk,:));

                        %% Set threshold
                        upper_th = ch_mean + 2*ch_std;
                        lower_th = ch_mean - 2*ch_std;

                        %% Find samples over or under amplitude threshold
                        I = obj.data.(tasks{ii,1}).LFP_info.region_LFP{jj,1}(kk,:) > upper_th | obj.data.(tasks{ii,1}).LFP_info.region_LFP{jj,1}(kk,:) < lower_th;
                        if sum(I)/length(I) > 0.2
                            obj.data.(tasks{ii,1}).LFP_info.bad_Chs{jj,1}(kk,1) = 1;
                        else
                            obj.data.(tasks{ii,1}).LFP_info.bad_Chs{jj,1}(kk,1) = 0;
                        end

                    end
                end
            end
        end

        %% Referencing
        function dataReference(obj,refScheme)

            %% function dataReference(obj,refScheme)

            %% Print
            fprintf('Running function dataReference..\n');

            %% Get the task names
            tasks = fieldnames(obj.data);

            %% Loop over tasks
            for ii = 1:length(tasks)

                %% Loop over wires
                for jj = 1:length(obj.data.(tasks{ii,1}).LFP_info.region_LFP)

                    %% Get the channel data and channel identification
                    ch_data = obj.data.(tasks{ii,1}).LFP_info.region_LFP{jj,1};
                    ch_bads = obj.data.(tasks{ii,1}).LFP_info.bad_Chs{jj,1};
                    ch_name = obj.data.(tasks{ii,1}).LFP_info.regionChan{jj,1};

                    %% Reject bad channels
                    ch_bads = logical(ch_bads);
                    ch_data = ch_data(~ch_bads,:);
                    ch_name = ch_name(~ch_bads,1);

                    %% Reference based on scheme
                    if strcmp(refScheme,'bipolar')

                        %% Loop to reference
                        for kk = length(ch_name):-1:1
                            if kk ~= 1
                                refData(kk,:) = ch_data(kk,:) - ch_data(kk-1,:);
                                refChan{kk,1} = append(ch_name{kk,1},'-',ch_name{kk-1,1});
                            end
                        end

                        %% Remove empty cells
                        I = cellfun('isempty',refChan);
                        refData = refData(~I,:);
                        refChan = refChan(~I,1);

                        %% Allocate to object
                        obj.data.(tasks{ii,1}).LFP_info.refData{jj,1} = refData;
                        obj.data.(tasks{ii,1}).LFP_info.refChan{jj,1} = refChan;

                        %% Clear
                        clear ch_bads ch_data ch_name refData refChan

                    elseif strcmp(refScheme,'laplace')

                        %% Loop to reference
                        for kk = length(ch_name):-1:1
                            if kk ~= length(ch_name) && kk ~= 1
                                refData(kk,:) = ch_data(kk,:) - ((ch_data(kk+1,:)+ch_data(kk-1,:)).*0.5);
                                refChan{kk,1} = append(ch_name{kk+1,1},'-',ch_name{kk,1},'-',ch_name{kk-1,1});
                            end
                        end

                        %% Remove empty cells
                        I = cellfun('isempty',refChan);
                        refData = refData(~I,:);
                        refChan = refChan(~I,1);

                        %% Allocate to object
                        obj.data.(tasks{ii,1}).LFP_info.refData{jj,1} = refData;
                        obj.data.(tasks{ii,1}).LFP_info.refChan{jj,1} = refChan;

                        %% Clear
                        clear ch_bads ch_data ch_name refData refChan

                    end
                end
            end
        end

        %% Filtering
        function dataFilter(obj,lfreq,hfreq,type,order)

            %% function dataFilter(obj,lfreq,hfreq,type,order)
            % This function filters the data in the object from lfreq to
            % hfreq using the filter type specified by 'type' and of user
            % specified order.  The filter type can be either 'fir' for a
            % standard FIR type filter, or 'iir' to implement a Butterworth
            % filter.

            %% Print
            fprintf('Running function dataFilter..\n');

            %% Get the task names
            tasks = fieldnames(obj.data);

            %% Loop over the tasks
            for ii = 1:length(tasks)

                %% Check for appropriate frequencies
                if exist('lfreq','var') && exist('hfreq','var')
                    if hfreq < lfreq
                        warning('Inappropriate hfreq.  Using sampling rate.\n')
                        hfreq = obj.data.(tasks{ii,1}).fs - 1;
                    end
                end

                if exist('lfreq','var') && exist('hfreq','var')
                    if lfreq > obj.data.(tasks{ii,1}).fs
                        warning('Inappropriate lfreq.  Setting to 1 Hz.\n')
                    end
                end

                if exist('lfreq','var') && ~exist('hfreq','var')
                    warning('No hfreq.  Using default value.\n')
                    hfreq = obj.data.(tasks{ii,1}).fs - 1;
                end

                if ~exist('lfreq','var') && exist('hfreq','var')
                    warning('No lfreq.  Using default value.')
                    lfreq = 1;
                end

                if ~exist('lfreq','var') && ~exist('hfreq','var')
                    warning('No lfreq or hfreq.  Using default values.')
                    lfreq = 1;
                    hfreq = obj.data.(tasks{ii,1}).fs - 1;
                end

                %% Loop over wires
                for jj = 1:length(obj.data.(tasks{ii,1}).LFP_info.refData)

                    %% Check filter type
                    if strcmp(type,'fir')

                        %% Create filter parameters
                        Wn = [lfreq,hfreq]./(obj.data.(tasks{ii,1}).fs./2);
                        b = fir1(order,Wn,'bandpass');

                        %% Filter
                        obj.data.(tasks{ii,1}).LFP_info.filtData{jj,1} = filtfilt(b,1,obj.data.(tasks{ii,1}).LFP_info.refData{jj,1}')';

                    elseif strcmp(type,'iir')

                        %% Create filter parameters
                        Wn = [lfreq,hfreq]./(obj.data.(tasks{ii,1}).fs./2);
                        [b,a] = butter(order,Wn,'bandpass');

                        %% Filter
                        obj.data.(tasks{ii,1}).LFP_info.filtData{jj,1} = filtfilt(b,a,obj.data.(tasks{ii,1}).LFP_info.refData{jj,1}')'; % Orig refData
                        
                    end
                end
            end
        end

        %% Identify TTL location in LFP
        function identifyTTLs(obj)

            %% function identifyTTL(obj)
            % This function identifies the LFP index corresponding to
            % pertinent events in the corresponding tasks and adds them to
            % the object.

            %% Print
            fprintf('Running function identifyTTL..\n');

            %% Get the task names
            tasks = fieldnames(obj.data);

            %% Loop over the tasks
            for ii = 1:length(tasks) %%%%%%%
                
                %% Find LFP indices for the events
                if strcmp(obj.task,'NEW_OLD_DELAY')
                    for jj = 1:length(obj.data.(tasks{ii,1}).info.NLX_TTLinfo.stimOn)
                        [~,obj.data.(tasks{ii,1}).eventIndex.stimOn(jj,1)] = min(abs(obj.data.(tasks{ii,1}).info.NLX_TTLinfo.stimOn(jj,1) - obj.data.(tasks{ii,1}).LFP_info.timestamps));
                        [~,obj.data.(tasks{ii,1}).eventIndex.stimOff(jj,1)] = min(abs(obj.data.(tasks{ii,1}).info.NLX_TTLinfo.stimOff(jj,1) - obj.data.(tasks{ii,1}).LFP_info.timestamps));
                        [~,obj.data.(tasks{ii,1}).eventIndex.question(jj,1)] = min(abs(obj.data.(tasks{ii,1}).info.NLX_TTLinfo.Prompt(jj,1) - obj.data.(tasks{ii,1}).LFP_info.timestamps));
                        [~,obj.data.(tasks{ii,1}).eventIndex.response(jj,1)] = min(abs(obj.data.(tasks{ii,1}).info.NLX_TTLinfo.responseTimes(jj,1) - obj.data.(tasks{ii,1}).LFP_info.timestamps));
                        [~,obj.data.(tasks{ii,1}).eventIndex.end(jj,1)] = min(abs(obj.data.(tasks{ii,1}).info.NLX_TTLinfo.endDelayPostResp(jj,1) - obj.data.(tasks{ii,1}).LFP_info.timestamps));
                    end

                    %% For baseline, find number of samples for 400 ms
                    nSamp = obj.data.(tasks{ii,1}).fs .* 0.400;
                    obj.data.(tasks{ii,1}).eventIndex.baselineStart = obj.data.(tasks{ii,1}).eventIndex.stimOn - nSamp;

                elseif strcmp(obj.task,'MEM_ORDER')

                    %% Grab the indices
                    obj.data.(tasks{ii,1}).eventIndex.crossOn = obj.data.(tasks{ii,1}).LFP_info.LFPIdx(obj.data.(tasks{ii,1}).LFP_info.eventIDs == 11);
                    obj.data.(tasks{ii,1}).eventIndex.clipOn = obj.data.(tasks{ii,1}).LFP_info.LFPIdx(obj.data.(tasks{ii,1}).LFP_info.eventIDs == 1);
                    obj.data.(tasks{ii,1}).eventIndex.clipOff = obj.data.(tasks{ii,1}).LFP_info.LFPIdx(obj.data.(tasks{ii,1}).LFP_info.eventIDs == 2);
                    obj.data.(tasks{ii,1}).eventIndex.question = obj.data.(tasks{ii,1}).LFP_info.LFPIdx(obj.data.(tasks{ii,1}).LFP_info.eventIDs == 3);

                    %% Create indices for crossOff and response
                    crossOff = zeros(length(obj.data.(tasks{ii,1}).eventIndex.crossOn),1);
                    response = zeros(length(obj.data.(tasks{ii,1}).eventIndex.crossOn),1);
                    for jj = 1:length(crossOff)
                        crossOff(jj,1) = obj.data.(tasks{ii,1}).behav.respMat(jj).CrossEnd - obj.data.(tasks{ii,1}).behav.respMat(jj).CrossStart;
                        response(jj,1) = obj.data.(tasks{ii,1}).behav.respMat(jj).respTime - obj.data.(tasks{ii,1}).behav.respMat(jj).QuesStart;
                    end

                    %% Modify from time and add to object
                    obj.data.(tasks{ii,1}).eventIndex.crossOff = obj.data.(tasks{ii,1}).eventIndex.crossOn + (floor(obj.data.(tasks{ii,1}).fs .* crossOff));
                    obj.data.(tasks{ii,1}).eventIndex.response = obj.data.(tasks{ii,1}).eventIndex.question + (floor(obj.data.(tasks{ii,1}).fs .* response));
                    
                end 
            end
        end

        %% Identify boundaries
        function identifyBoundaries(obj)

            %% function identifyBoundaries(obj)
            % This function grabs the clip boundaries for the MemOrder
            % task.

            %% Check if subject needs boundaries
            if strcmp(obj.task,'MEM_ORDER')

                %% Print
                fprintf('Running function identifyBoundaries..\n')

                %% Get task names
                tasks = fieldnames(obj.data);

                %% Loop over tasks
                for ii = 1:length(tasks)
                    if strcmp(tasks{ii,1},'Encode')
                        obj.data.(tasks{ii,1}).boundaries = extractBefore({obj.data.(tasks{ii,1}).behav.respMat.ClipName}','_');
                    else
                        obj.data.(tasks{ii,1}).boundaries = extractBefore({obj.data.(tasks{ii,1}).behav.respMat.FrameName}','_');
                    end
                end

            end
        end

        %% Allocate data
        function allocateData(obj)

            %% function allocateData(obj)
            % This function chucks the data into appropriate segments from
            % the LFP based on the task that was performed.

            %% Print
            fprintf('Running function allocateData..\n');

            %% Get task names
            tasks = fieldnames(obj.data);

            %% Loop over tasks
            for ii = 1:length(tasks)

                %% Remove trials that are NaN
                if strcmp(obj.task,'NEW_OLD_DELAY')
                    I = isnan(obj.data.(tasks{ii,1}).info.NLX_TTLinfo.stimOn);
                    obj.data.(tasks{ii,1}).eventIndex.stimOn = obj.data.(tasks{ii,1}).eventIndex.stimOn(~I);
                    obj.data.(tasks{ii,1}).eventIndex.stimOff = obj.data.(tasks{ii,1}).eventIndex.stimOff(~I);
                    obj.data.(tasks{ii,1}).eventIndex.question = obj.data.(tasks{ii,1}).eventIndex.question(~I);
                    obj.data.(tasks{ii,1}).eventIndex.response = obj.data.(tasks{ii,1}).eventIndex.response(~I);
                    obj.data.(tasks{ii,1}).eventIndex.baselineStart = obj.data.(tasks{ii,1}).eventIndex.baselineStart(~I);
                    obj.data.(tasks{ii,1}).eventIndex.end = obj.data.(tasks{ii,1}).eventIndex.end(~I);
                    clear I
                end

                %% Loop over wires
                for jj = 1:length(obj.data.(tasks{ii,1}).LFP_info.filtData)

                    %% Write anonymous function
                    extractSegments = @(s,e) obj.data.(tasks{ii,1}).LFP_info.filtData{jj,1}(:,s:e-1);

                    %% Extracy segments
                    if strcmp(obj.task,'NEW_OLD_DELAY')
                        obj.data.(tasks{ii,1}).eventData.baseline{jj,1} = arrayfun(@(i) extractSegments(obj.data.(tasks{ii,1}).eventIndex.baselineStart(i), obj.data.(tasks{ii,1}).eventIndex.stimOn(i)), 1:length(obj.data.(tasks{ii,1}).eventIndex.stimOn), 'UniformOutput', false)';
                        obj.data.(tasks{ii,1}).eventData.presentation{jj,1} = arrayfun(@(i) extractSegments(obj.data.(tasks{ii,1}).eventIndex.stimOn(i), obj.data.(tasks{ii,1}).eventIndex.stimOff(i)), 1:length(obj.data.(tasks{ii,1}).eventIndex.stimOn), 'UniformOutput', false)';
                        obj.data.(tasks{ii,1}).eventData.question{jj,1} = arrayfun(@(i) extractSegments(obj.data.(tasks{ii,1}).eventIndex.question(i), obj.data.(tasks{ii,1}).eventIndex.response(i)), 1:length(obj.data.(tasks{ii,1}).eventIndex.stimOn), 'UniformOutput', false)';
                        obj.data.(tasks{ii,1}).eventData.epoch{jj,1} = arrayfun(@(i) extractSegments(obj.data.(tasks{ii,1}).eventIndex.baselineStart(i), obj.data.(tasks{ii,1}).eventIndex.response(i)), 1:length(obj.data.(tasks{ii,1}).eventIndex.stimOn), 'UniformOutput', false)';
                    else
                        obj.data.(tasks{ii,1}).eventData.baseline{jj,1} = arrayfun(@(i) extractSegments(obj.data.(tasks{ii,1}).eventIndex.crossOn(i), obj.data.(tasks{ii,1}).eventIndex.crossOff(i)), 1:length(obj.data.(tasks{ii,1}).eventIndex.crossOn), 'UniformOutput', false)';
                        obj.data.(tasks{ii,1}).eventData.presentation{jj,1} = arrayfun(@(i) extractSegments(obj.data.(tasks{ii,1}).eventIndex.clipOn(i), obj.data.(tasks{ii,1}).eventIndex.clipOff(i)), 1:length(obj.data.(tasks{ii,1}).eventIndex.crossOn), 'UniformOutput', false)';
                        obj.data.(tasks{ii,1}).eventData.question{jj,1} = arrayfun(@(i) extractSegments(obj.data.(tasks{ii,1}).eventIndex.question(i), obj.data.(tasks{ii,1}).eventIndex.response(i)), 1:length(obj.data.(tasks{ii,1}).eventIndex.crossOn), 'UniformOutput', false)';
                        obj.data.(tasks{ii,1}).eventData.epoch{jj,1} = arrayfun(@(i) extractSegments(obj.data.(tasks{ii,1}).eventIndex.crossOn(i), obj.data.(tasks{ii,1}).eventIndex.response(i)), 1:length(obj.data.(tasks{ii,1}).eventIndex.crossOn), 'UniformOutput', false)';
                    end
                end
            end
        end

        %% Baseline correct and Z-score
        function baselineCorrectZscore(obj)

            %% function baselineCorrectZscore(obj)
            % This function performs baseline correction using the baseline
            % period, and z-score normalized the LFP.

            %% Print
            fprintf('Running function baselineCorrectZscore..\n')

            %% Get the tasks
            tasks = fieldnames(obj.data);

            %% Loop over tasks
            for ii = 1:length(tasks)

                %% Loop over wire
                for jj = 1:length(obj.data.(tasks{ii,1}).eventData.baseline)

                    %% Loop over epochs
                    for kk = 1:length(obj.data.(tasks{ii,1}).eventData.baseline{jj,1})

                        %% Calculate channel means from baseline
                        ch_mean = mean(obj.data.(tasks{ii,1}).eventData.baseline{jj,1}{kk,1},2);

                        %% Subtract channel means
                        obj.data.(tasks{ii,1}).eventData.baseline{jj,1}{kk,1} = obj.data.(tasks{ii,1}).eventData.baseline{jj,1}{kk,1} - ch_mean;
                        obj.data.(tasks{ii,1}).eventData.presentation{jj,1}{kk,1} = obj.data.(tasks{ii,1}).eventData.presentation{jj,1}{kk,1} - ch_mean;
                        obj.data.(tasks{ii,1}).eventData.question{jj,1}{kk,1} = obj.data.(tasks{ii,1}).eventData.question{jj,1}{kk,1} - ch_mean;
                        obj.data.(tasks{ii,1}).eventData.epoch{jj,1}{kk,1} = obj.data.(tasks{ii,1}).eventData.epoch{jj,1}{kk,1} - ch_mean;

                        %% Z-score the data
                        obj.data.(tasks{ii,1}).eventData.baseline{jj,1}{kk,1} = normalize(obj.data.(tasks{ii,1}).eventData.baseline{jj,1}{kk,1},2,"zscore");
                        obj.data.(tasks{ii,1}).eventData.presentation{jj,1}{kk,1} = normalize(obj.data.(tasks{ii,1}).eventData.presentation{jj,1}{kk,1},2,"zscore");
                        obj.data.(tasks{ii,1}).eventData.question{jj,1}{kk,1} = normalize(obj.data.(tasks{ii,1}).eventData.question{jj,1}{kk,1},2,"zscore");
                        obj.data.(tasks{ii,1}).eventData.epoch{jj,1}{kk,1} = normalize(obj.data.(tasks{ii,1}).eventData.epoch{jj,1}{kk,1},2,"zscore");
                        
                        %% Clear
                        clear ch_mean;

                    end
                end
            end
        end

        %% Plot phases
        function plotPhases(obj)

            %% function plotPhases(obj)
            % This function acts as a sanity check to make sure the
            % indices for the task don't overlap, and to give an overview
            % of what the LFP data look like.

            %% Print
            fprintf('Running function plotPhases..\n')

            %% Get the tasks
            tasks = fieldnames(obj.data);

            %% Loop over tasks
            for ii = 1:length(tasks)

                %% Create title
                if strcmp(obj.task,'NEW_OLD_DELAY')
                    title_str = append(extractBefore(obj.subID,'_'),' ','Task',' ',tasks{ii,1});
                else
                    title_str = append(obj.subID,' ','Task',' ',tasks{ii,1});
                end

                %% Get data
                filtData = cell2mat(obj.data.(tasks{ii,1}).LFP_info.filtData);

                %% Plot LFP data
                figure(ii);
                fig = gca;
                fig.XAxis.FontSize = 24;
                fig.XAxis.FontWeight = 'bold';
                fig.YAxis.FontSize = 24;
                fig.YAxis.FontWeight = 'bold';
                fig.XLim = ([min(obj.data.(tasks{ii,1}).LFP_info.converted_timestamps) max(obj.data.(tasks{ii,1}).LFP_info.converted_timestamps)]);
                fig.YLim = ([1.1*(min(min(filtData))) 1.1*max(max(filtData))]);
                title(title_str,'FontWeight','bold','FontSize',28);
                hold on
                for z = 1:size(filtData,1)
                    plot(obj.data.(tasks{ii,1}).LFP_info.converted_timestamps,filtData(z,:))
                end

                %% Figure properties
                xlabel('Time (sec)','FontWeight','bold','FontSize',24);
                ylabel('Amplitude (V)','FontWeight','bold','FontSize',24);

                %% Bands for baseline
                if strcmp(obj.task,'NEW_OLD_DELAY')
                    bands = [obj.data.(tasks{ii,1}).eventIndex.baselineStart,obj.data.(tasks{ii,1}).eventIndex.stimOn]./obj.data.(tasks{ii,1}).fs;
                else
                    bands = [obj.data.(tasks{ii,1}).eventIndex.crossOn + obj.data.(tasks{ii,1}).LFP_info.correction_factor,obj.data.(tasks{ii,1}).eventIndex.crossOff+ obj.data.(tasks{ii,1}).LFP_info.correction_factor]./obj.data.(tasks{ii,1}).fs;
                end
                xp = [bands fliplr(bands)];
                yp = ([[1;1]*1.1*min(ylim); [1;1]*1.1*max(ylim)]*ones(1,size(bands,1))).';
                for k = 1:size(bands,1)
                    patch(xp(k,:),yp(k,:),[1 0 0],'FaceAlpha',0.1,'EdgeColor','none')
                end

                %% Bands for presentation
                if strcmp(obj.task,'NEW_OLD_DELAY')
                    bands = [obj.data.(tasks{ii,1}).eventIndex.stimOn,obj.data.(tasks{ii,1}).eventIndex.stimOff]./obj.data.(tasks{ii,1}).fs;
                else
                    bands = [obj.data.(tasks{ii,1}).eventIndex.clipOn + obj.data.(tasks{ii,1}).LFP_info.correction_factor,obj.data.(tasks{ii,1}).eventIndex.clipOff + obj.data.(tasks{ii,1}).LFP_info.correction_factor]./obj.data.(tasks{ii,1}).fs;
                end
                xp = [bands fliplr(bands)];
                yp = ([[1;1]*1.1*min(ylim); [1;1]*1.1*max(ylim)]*ones(1,size(bands,1))).';
                for k = 1:size(bands,1)
                    patch(xp(k,:),yp(k,:),[0 1 0],'FaceAlpha',0.1,'EdgeColor','none')
                end

                %% Bands for question
                if strcmp(obj.task,'NEW_OLD_DELAY')
                    bands = [obj.data.(tasks{ii,1}).eventIndex.question,obj.data.(tasks{ii,1}).eventIndex.response]./obj.data.(tasks{ii,1}).fs;
                else
                    bands = [obj.data.(tasks{ii,1}).eventIndex.question + obj.data.(tasks{ii,1}).LFP_info.correction_factor,obj.data.(tasks{ii,1}).eventIndex.response + obj.data.(tasks{ii,1}).LFP_info.correction_factor]./obj.data.(tasks{ii,1}).fs;
                end
                xp = [bands fliplr(bands)];
                yp = ([[1;1]*1.1*min(ylim); [1;1]*1.1*max(ylim)]*ones(1,size(bands,1))).';
                for k = 1:size(bands,1)
                    patch(xp(k,:),yp(k,:),[0 0 1],'FaceAlpha',0.1,'EdgeColor','none')
                end

            end
        end

        %% Saccade analysis
        function saccadeAnalysis(obj,lfreq,hfreq,type,order)

            %% function saccadeAnalysis(obj,lfreq,hfreq,type,order)
            % This function examines the saccades in filtered data and
            % calculated the inter-trial phase coherence to examine whether
            % or not saccades resulted in a phase reset.

            %% Print
            fprintf('Running function saccadeAnalysis..\n');

            %% Get the tasks
            tasks = fieldnames(obj.data);
            
            %% Loop over the tasks
            for ii = 1:length(tasks)

                %% Create filter parameters
                Wn = [lfreq,hfreq]./(obj.data.(tasks{ii,1}).fs/2);
                if strcmp(type,'fir')
                    b = fir1(order,Wn,'bandpass');
                    a = 1;
                else
                    [b,a] = butter(order,Wn,'bandpass');
                end

                %% Set nIterations for statistical testing
                nIter = 1000; % 1000 for paper

                %% Calculate the number of sample pre and post saccade
                nSamp = (obj.data.(tasks{ii,1}).fs*400)/1000; % 400 ms pre and post

                %% Get the eye and ttl info
                if strcmp(obj.task,'NEW_OLD_DELAY')
                    GAZEcl = obj.data.(tasks{ii,1}).info.GazeEVENTs{:,end};
                    TTLinfo = obj.data.(tasks{ii,1}).info.GazeEVENTs{:,19};
                else
                    eye = obj.data.(tasks{ii,1}).eye2use;
                    fields = fieldnames(obj.data.(tasks{ii,1}).eyeInfo);
                    I = contains(fields,eye,'IgnoreCase',true);
                    eyeField = fields(I);

                    GAZEcl = obj.data.(tasks{ii,1}).eyeInfo.(eyeField{1,1}){:,end};
                    TTLinfo = obj.data.(tasks{ii,1}).eyeInfo.TTLinfo;
                end

                %% Set parameters
                nEpochs = length(GAZEcl);
                nWire = length(obj.data.(tasks{ii,1}).LFP_info.filtData);

                %% Loop over wires
                for jj = 1:nWire
                    jj

                    %% Grab the data
                    tmpData = obj.data.(tasks{ii,1}).LFP_info.filtData{jj,1}';

                    %% Get the number of channels on the wire
                    nChan = size(tmpData,2);

                    %% Filter the data
                    tmpData = filtfilt(b,a,tmpData);

                    %% Calculate phase
                    phase = angle(hilbert(tmpData));

                    %% Loop over electrodes
                    for kk = 1:nChan

                        %% Calculate CWT
                        %[p,f] = cwt(tmpData(:,kk),obj.data.(tasks{ii,1}).fs);

                        %% Preallocate
                        randSaccade = cell(nEpochs,1);

                        %% Loop over epochs
                        for ll = 1:nEpochs

                            %% Grab the NLX and ELNK times
                            NLXstart = double(TTLinfo{ll,1}{2,4});
                            ELNKstart = double(TTLinfo{ll,1}{2,5});

                            %% Check if saccades occurred
                            if istable(GAZEcl{ll,1}.saccades)

                                %% Grab the saccade times
                                saccades = GAZEcl{ll,1}.saccades;
                                saccadeOnset = double(saccades.starttime);

                                %% Saccade duration
                                duration = double(saccades.duration);

                                %% Calculate NLX sample difference
                                saccadeOnset = saccadeOnset - ELNKstart;
                                saccadeOnset = floor((saccadeOnset./1000) .* obj.data.(tasks{ii,1}).fs);

                                %% Identify NLX onset, start, and stop indices
                                if strcmp(obj.task,'MEM_ORDER')
                                    saccadeOnset = saccadeOnset - obj.data.(tasks{ii,1}).LFP_info.correction_factor + NLXstart;
                                    saccadeStart = saccadeOnset - nSamp;
                                    saccadeEnd = saccadeOnset + nSamp - 1;
                                else
                                    [~,I] = min(abs(NLXstart - obj.data.(tasks{ii,1}).LFP_info.timestamps));
                                    saccadeOnset = saccadeOnset + I;
                                    saccadeStart = saccadeOnset - nSamp;
                                    saccadeEnd = saccadeOnset + nSamp - 1;
                                end

                                %% Loop over saccades
                                for mm = 1:length(saccadeStart)

                                    %% Skip saccade if it goes beyond the recording
                                    if saccadeStart(mm,1) > 0 && saccadeEnd(mm,1) <= length(phase)

                                        %% Grab the saccade phases
                                        tmpSaccade = phase(saccadeStart(mm,1):saccadeEnd(mm,1),kk);
                                        saccadePhase{jj,1}{kk,1}{ll,1}(mm,1:saccadeEnd(mm,1)-saccadeStart(mm,1)+1) = tmpSaccade;

                                        %% Grab the saccade duration
                                        saccadeDuration{jj,1}{kk,1}{ll,1}(mm,1) = duration(mm,1);

                                        %% Loop for statistical testing
                                        for nn = 1:nIter

                                            %% Randomly permute
                                            ix = randperm(length(tmpSaccade));
                                            %randSaccade{jj,1}{kk,1}{ll,1}(mm,:,nn) = tmpSaccade(ix);
                                            randSaccade{ll,1}(mm,:,nn) = tmpSaccade(ix);

                                        end

                                        %% Clear
                                        clear tmpSaccade
                                    else
                                        saccadePhase{jj,1}{kk,1}{ll,1}(mm,1:saccadeEnd(mm,1)-saccadeStart(mm,1)+1) = NaN;
                                    end

                                end % End saccades

                                %% Clear
                                clear saccades saccadeStart

                            end

                            %% Clear
                            clear NLXstart ELNKstart
                        end % End epochs

                        %% Remove empty channels
                        randSaccade = randSaccade(~cellfun('isempty',randSaccade));

                        %% Calculate ITPC
                        for iter = 1:nIter
                            for tmpEpoch = 1:length(randSaccade)
                                sacITPC{tmpEpoch,1} = randSaccade{tmpEpoch,1}(:,:,iter);
                            end
                            sacITPC = cell2mat(sacITPC);
                            for z = 1:size(sacITPC,2)
                                iterITPC{jj,1}{kk,1}(iter,z) = abs(mean(exp(1i.*sacITPC(:,z))));
                            end
                            clear sacITPC
                        end

                        %% Clear
                        clear p f randSaccade;

                    end % End channels
                end % End wires

                %% Can loop here to calculate means if desired

                %% Loop to plot saccades
                for jj = 1:nWire
                    for kk = 1:length(saccadePhase{jj,1})

                        %% Get the data
                        saccadeData = cell2mat(saccadePhase{jj,1}{kk,1});
                        ch_name = obj.data.(tasks{ii,1}).LFP_info.refChan{jj,1}{kk,1};
                        saccade_sample = -200:1:199;
                        saccade_number = 1:1:size(saccadeData,1);

                        %% Create the figure
                        figure('Position',[1,49,2560,1315])
                        tl = tiledlayout(10,1);
                        txt = title(tl,ch_name);
                        txt.FontWeight = 'bold';
                        txt.FontSize = 32;
                        xlabel(tl,'Sample Number','FontSize',20,'FontWeight','bold');

                        % Tile 1
                        nt = nexttile([7 1]);
                        I = any(isnan(saccadeData),2);
                        saccadeData = saccadeData(~I,:);
                        imagesc(saccade_sample,saccade_number,saccadeData);
                        xticks([]);
                        ylabel('Saccade Number','FontSize',20,'FontWeight','bold')
                        nt.FontWeight = 'bold';
                        nt.FontSize = 20;

                        % Calculate ITPC
                        itpc = NaN(1,size(saccadeData,2));
                        for z = 1:length(itpc)
                            itpc(1,z) = abs(mean(exp(1i.*saccadeData(:,z))));
                        end

                        % Grab random phases
                        randPhase = iterITPC{jj,1}{kk,1};

                        % Calculate value for 95-percentile
                        P = prctile(randPhase,95,1);

                        % Find where ITPC is greater than random
                        I = itpc > P;

                        % Get the bands
                        samps = saccade_sample(1,I);
                        samp_start = samps - 0.5;
                        samp_end = samps + 0.5;
                        bands = [samp_start',samp_end'];
                        

                        % Tile 2
                        nt = nexttile([3 1]);
                        plot(saccade_sample,itpc,'LineWidth',3,'Color','black')
                        xp = [bands fliplr(bands)];
                        %yp = ([[0;0]*1.1*min(ylim); [1;1]*1.1*max(ylim)]*ones(1,size(bands,1))).';
                        yp = ([[1;1]*1.1*min(ylim); [1;1]*1.1*max(ylim)]*ones(1,size(bands,1))).';
                        for k = 1:size(bands,1)
                            patch(xp(k,:),yp(k,:),[0.7 0.7 0.7],'FaceAlpha',0.1,'EdgeColor','none')
                        end
                        xlim([-200 200]);
                        ylim([0 1.1*max(itpc)]);
                        ylabel('ITPC','FontSize',20,'FontWeight','bold')
                        nt.FontWeight = 'bold';
                        nt.FontSize = 20;

                        % clear
                        clear saccadeData itpc samps samp_start samp_end randPhase P I

                    end
                end
                x = 0;
            end
            x = 0;

        end

        %% Next function


    end
end