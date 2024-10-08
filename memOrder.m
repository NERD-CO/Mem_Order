classdef memOrder < handle

    %% Properties
    properties(SetAccess = private)
        subID
        folderDir
        encode
        sceneRecog
        timeDiscrim
        validTasks
        validPhases
    end

    %% Methods
    methods

        %% Constructor
        function obj = memOrder(folderDir)

            %% Notes

            %% Set folderDir and subID
            obj.subID = extractBefore(extractAfter(folderDir,'Folder\'),'\');
            obj.folderDir = folderDir;

            %% Set the directories
            behavDir = append(folderDir,'Behavioral_Data\Raw\');
            eyeDir = append(folderDir,'Eye-tracking\Raw\');
            nwbDir = append(folderDir,'NWBProcessing\NWB_Data\');

            %% Determine file IDs
            info_file = dir(fullfile(folderDir,'*.xlsx'));
            info_file = append(folderDir,info_file.name);
            fileIDs = readtable(info_file);

            %% Grab file names
            nwb_names = fileIDs.File;
            behav_names = fileIDs.Behav;
            eye_names = fileIDs.Eye;

            %% Create logicals
            encode = logical(fileIDs.Encode);
            sceneRecog = logical(fileIDs.SceneRecognition);
            timeDiscrim = logical(fileIDs.TimeDiscrimination);

            %% Add to object if task was performed
            if any(encode == 1)
                obj.encode.nwb_fname = append(nwbDir,nwb_names{encode});
                obj.encode.behave_fname = append(behavDir,behav_names{encode});
                obj.encode.eye_fname = append(eyeDir,eye_names{encode});
            end

            if any(sceneRecog == 1)
                obj.sceneRecog.nwb_fname = append(nwbDir,nwb_names{sceneRecog});
                obj.sceneRecog.behave_fname = append(behavDir,behav_names{sceneRecog});
                obj.sceneRecog.eye_fname = append(eyeDir,eye_names{sceneRecog});
            end

            if any(timeDiscrim == 1)
                obj.timeDiscrim.nwb_fname = append(nwbDir,nwb_names{timeDiscrim});
                obj.timeDiscrim.behave_fname = append(behavDir,behav_names{timeDiscrim});
                obj.timeDiscrim.eye_fname = append(eyeDir,eye_names{timeDiscrim});
            end

            %% Find the non-empty structures for dynamic indexing
            fieldNames = fieldnames(obj);
            matchedNames = cell(length(fieldNames),1);
            for ii = 1:length(fieldNames)
                if isstruct(obj.(fieldNames{ii,1}))
                    if ~isempty(obj.(fieldNames{ii,1}))
                        matchedNames{ii,1} = fieldNames{ii,1};
                    end
                end
            end

            %% Remove empty matchedNames
            matchedNames = matchedNames(~cellfun('isempty',matchedNames));

            %% Add to object
            obj.validTasks = matchedNames;

            %% Preallocate
            tmpNames = cell(length(matchedNames),1);
            count = 1;

            %% Loop through and do the NWB processing
            for ii = 1:length(matchedNames)

                %% Compare names to identify times through the same file
                if any(strcmp(obj.(matchedNames{ii,1}).nwb_fname,tmpNames))
                    count = count + 1;
                else
                    count = 1;
                end

                %% Move the file name to tmpNames
                tmpNames{ii,1} = obj.(matchedNames{ii,1}).nwb_fname;

                %% Print
                fprintf('Ephys data loaded for sub %s and task %s..\n',obj.subID,matchedNames{ii,1})

                %% Load the data
                tmp = nwbRead(obj.(matchedNames{ii}).nwb_fname);

                %% Load the behavioral info here
                obj.(matchedNames{ii}).behave_info = load(obj.(matchedNames{ii}).behave_fname);

                %% Convert respmat vales from old (-3 to 3) to new (111 to 116)

                %% Load macrowire time stamps
                timestamps = tmp.processing.get('ecephys').nwbdatainterface.get...
                    ('LFP').electricalseries.get('MacroWireSeries').timestamps.load;

                %% Downsample time stamps
                obj.(matchedNames{ii,1}).LFP_timestamps = downsample(timestamps,8);

                %% Convert time stamps
                obj.(matchedNames{ii,1}).LFP_converted_timestamps = (obj.(matchedNames{ii,1}).LFP_timestamps - obj.(matchedNames{ii,1}).LFP_timestamps(1,1))/(1e6);

                %% Get the sampling frequency
                LFP_sessionInfo = tmp.processing.get('ecephys').nwbdatainterface.get...
                    ('LFP').electricalseries.get('MacroWireSeries');
                obj.(matchedNames{ii,1}).fs = str2double(cell2mat(extractBetween(LFP_sessionInfo.description,'= ',':')));

                %% Get voltage for macrowires and convert
                obj.(matchedNames{ii,1}).LFP_data = (double(tmp.processing.get('ecephys').nwbdatainterface.get...
                    ('LFP').electricalseries.get('MacroWireSeries').data.load)) .* LFP_sessionInfo.data_conversion;

                %% Channel IDs
                chanLabels = cellstr(tmp.general_extracellular_ephys_electrodes.vectordata.get('label').data.load()); %use MA only!
                MAchan = find(contains(chanLabels,'MA_'));
                chanID = cellstr(tmp.general_extracellular_ephys_electrodes.vectordata.get('location').data.load());
                hemisphere = cellstr(tmp.general_extracellular_ephys_electrodes.vectordata.get('hemisph').data.load());
                shortBnames = cellstr(tmp.general_extracellular_ephys_electrodes.vectordata.get('shortBAn').data.load());
                wireID = tmp.general_extracellular_ephys_electrodes.vectordata.get('channID').data.load();

                obj.(matchedNames{ii,1}).chanID = chanID(MAchan);
                obj.(matchedNames{ii,1}).chanHemi = hemisphere(MAchan);
                obj.(matchedNames{ii,1}).chanSname = shortBnames(MAchan);
                obj.(matchedNames{ii,1}).wireID = wireID(MAchan);

                %% Find unique brain regions and turn into a table
                tmpBRegUni = unique(obj.(matchedNames{ii,1}).chanSname);
                hemiTemp = cell(length(tmpBRegUni),1);
                longName = cell(length(tmpBRegUni),1);
                for ui = 1:length(tmpBRegUni)
                    tmpIND = find(matches(obj.(matchedNames{ii,1}).chanSname,tmpBRegUni{ui}),1,'first');
                    hemiTemp{ui} = obj.(matchedNames{ii,1}).chanHemi{tmpIND};
                    longName{ui} = obj.(matchedNames{ii,1}).chanID{tmpIND};
                end
                obj.(matchedNames{ii,1}).brTABLE = table(tmpBRegUni,hemiTemp,longName,'VariableNames',{'SEEGele',...
                    'Hemisphere','LongBRname'});

                %% Assign LFP to each brain region
                obj.(matchedNames{ii,1}).regionLFP = cell(length(obj.(matchedNames{ii,1}).brTABLE.SEEGele),1);
                for k = 1:height(obj.(matchedNames{ii,1}).brTABLE)
                    reg = strcmp(obj.(matchedNames{ii,1}).brTABLE.SEEGele{k},obj.(matchedNames{ii,1}).chanSname);
                    obj.(matchedNames{ii,1}).regionLFP{k} = obj.(matchedNames{ii,1}).LFP_data(reg,:);
                end

                %% Extract event key
                obj.(matchedNames{ii,1}).eventTimes = tmp.acquisition.get('events').timestamps.load();
                temp_eventIDs = cellstr(tmp.acquisition.get('events').data.load());

                %% Convert eventIDs from hexadecimal
                I = contains(temp_eventIDs,'TTL');
                obj.(matchedNames{ii,1}).eventTimes = obj.(matchedNames{ii,1}).eventTimes(I);
                temp_eventIDs2 = temp_eventIDs(I);
                TTL_task = extractBetween(temp_eventIDs2,'(',')');
                obj.(matchedNames{ii,1}).eventIDs = cellfun(@(x) hex2dec(x),TTL_task,'UniformOutput',true);

                %% Find task start and stop indices
                x = find(obj.(matchedNames{ii,1}).eventIDs == 61); % task start
                y = find(obj.(matchedNames{ii,1}).eventIDs == 60); % task end

                %% Find the time in the recording that corresponds to the indices
                for jj = 1:length(obj.(matchedNames{ii,1}).eventTimes)
                    [~,obj.(matchedNames{ii,1}).LFPIdx(jj,1)] = min(abs(obj.(matchedNames{ii,1}).eventTimes(jj,1) - obj.(matchedNames{ii,1}).LFP_timestamps));
                end

                %% Little fix for uneven indices if needed
                if length(y) < length(x) % if no stop
                    a = length(x);
                    y(a,1) = length(obj.(matchedNames{ii,1}).eventIDs);
                end

                if length(x) < length(y) % if no start
                    x(1,1) = 1;
                end

                %% Reduce indices to those relevant to the task
                obj.(matchedNames{ii,1}).eventIDs = obj.(matchedNames{ii,1}).eventIDs(x(count,1):y(count,1),1);
                obj.(matchedNames{ii,1}).eventTimes = obj.(matchedNames{ii,1}).eventTimes(x(count,1):y(count,1),1);
                obj.(matchedNames{ii,1}).LFPIdx = obj.(matchedNames{ii,1}).LFPIdx(x(count,1):y(count,1),1);

                %% Reduce data to that of the individual task
                obj.(matchedNames{ii,1}).LFP_timestamps = obj.(matchedNames{ii,1}).LFP_timestamps(obj.(matchedNames{ii,1}).LFPIdx(1,1):obj.(matchedNames{ii,1}).LFPIdx(end,1),:);
                obj.(matchedNames{ii,1}).LFP_converted_timestamps = obj.(matchedNames{ii,1}).LFP_converted_timestamps(obj.(matchedNames{ii,1}).LFPIdx(1,1):obj.(matchedNames{ii,1}).LFPIdx(end,1),:);
                obj.(matchedNames{ii,1}).LFP_data = obj.(matchedNames{ii,1}).LFP_data(:,obj.(matchedNames{ii,1}).LFPIdx(1,1):obj.(matchedNames{ii,1}).LFPIdx(end,1));

                %% Calculate correction factor
                obj.(matchedNames{ii,1}).correction_factor = obj.(matchedNames{ii,1}).LFPIdx(1,1) - 1;
                obj.(matchedNames{ii,1}).LFPIdx = obj.(matchedNames{ii,1}).LFPIdx - obj.(matchedNames{ii,1}).correction_factor;

                %% Find respMat table
                % fieldNames = fieldnames(obj.(matchedNames{ii}).respMat);
                % respMat_name = cell(length(fieldNames),1);
                % for z = 1:length(fieldNames)
                %     if contains(fieldNames{z,1},'respMat')
                %         respMat_name{z,1} = fieldNames{z,1};
                %     end
                % end

                %% Remove empty cells
                % respMat_name = respMat_name(~cellfun('isempty',respMat_name));

                %% Convert response TTLs if they exist
                % responseVals = [111 112 113 114 115 116];
                % originalVals = [-3 -2 -1 1 2 3];
                % for z = 1:length(obj.(matchedNames{ii,1}).respMat.(respMat_name{1,1}))
                %     [~,idx] = ismember(obj.(matchedNames{ii,1}).respMat.(respMat_name{1,1})(z).respValue,originalVals);
                %     if ~isempty(idx)
                %         obj.(matchedNames{ii,1}).confidenceRating(z,1) = responseVals(idx);
                %     else
                %         obj.(matchedNames{ii,1}).confidenceRating(z,1) = NaN;
                %     end
                % end

                %% Test
                % count1 = 1;
                % for a = 1:length(obj.(matchedNames{ii,1}).eventIDs)
                %     if (obj.(matchedNames{ii,1}).eventIDs(a,1) < 111) || (isnan(obj.(matchedNames{ii,1}).eventIDs(a,1)))
                %         continue
                %     elseif (obj.(matchedNames{ii,1}).eventIDs(a,1) >= 111) && (obj.(matchedNames{ii,1}).eventIDs(a+1,1) >= 111)
                %         num1 = obj.(matchedNames{ii,1}).eventIDs(a,1);
                %         num2 = obj.(matchedNames{ii,1}).eventIDs(a+1,1);
                %         confRating = obj.(matchedNames{ii,1}).confidenceRating(count1,1);
                %         if num1 == confRating
                %             obj.(matchedNames{ii,1}).eventIDs(a+1,1) = NaN;
                %             count1 = count1 + 1;
                %         elseif num2 == confRating
                %             obj.(matchedNames{ii,1}).eventIDs(a,1) = NaN;
                %             count1 = count1 + 1;
                %         end
                %     end
                % end

                %% Find NaNs if they exist
                % I = isnan(obj.(matchedNames{ii,1}).eventIDs);

            end
        end

        %% Select channels
        function selectChannels(obj,region,hemi)

            %% Notes
            % This function loops through valid tasks to find the channels
            % of interest for the analysis.

            %% Print
            fprintf('Running function selectChannels..\n');

            %% Loop through tasks and select channels
            for ii = 1:length(obj.validTasks)

                %% Get the task
                task = obj.validTasks{ii,1};

                %% Find channels of interest
                I = strcmp(obj.(task).chanID,region) & strcmp(obj.(task).chanHemi,hemi);

                %% Remove channels
                obj.(task).LFP_data = obj.(task).LFP_data(I,:);
                obj.(task).chanID = obj.(task).chanID(I,1);
                obj.(task).chanHemi = obj.(task).chanHemi(I,1);
                obj.(task).chanSname = obj.(task).chanSname(I,1);
                obj.(task).wireID = obj.(task).wireID(I,1);

                %% Modify channel names
                for z = 1:length(obj.(task).chanID)
                    obj.(task).chanIdx{z,1} = strcat(obj.(task).chanSname{z,1},'_',num2str(z));
                end
            end
        end

        %% Identify bad channels
        function identifyBads(obj)

            %% Notes
            % This function loops through valid tasks and find individual
            % channels in each task that exceed a certain threshold.  Those
            % that do are replaced with NaN's and are not used in
            % subsequent functions.

            %% Print
            fprintf('Running function identifyBads..\n');

            %% Loop through tasks
            for ii = 1:length(obj.validTasks)

                %% Get the task
                task = obj.validTasks{ii,1};

                %% Loop through channels to identify bads
                for jj = 1:length(obj.(task).chanIdx)

                    %% Calculate mean and standard deviation
                    ch_mean = mean(obj.(task).LFP_data(jj,:));
                    ch_std = std(obj.(task).LFP_data(jj,:));

                    %% Set thresholds
                    upper_th = ch_mean + 2*ch_std;
                    lower_th = ch_mean - 2*ch_std;

                    %% Determine if a channel is bad by comparing to its thresholds
                    if sum(obj.(task).LFP_data(jj,:) > upper_th | obj.(task).LFP_data(jj,:) < lower_th) > 0.2*length(obj.(task).LFP_data(jj,:))
                        obj.(task).keepIdx(jj,1) = 0;
                    else
                        obj.(task).keepIdx(jj,1) = 1;
                    end
                end

                %% Convert to logical
                obj.(task).keepIdx = logical(obj.(task).keepIdx);

                %% Set bad channels to NaN
                obj.(task).LFP_data(~obj.(task).keepIdx,:) = NaN;
            end
        end

        %% High pass filter
        function dataFilter(obj,lfreq)

            %% Notes
            % This function high-passes the LFP data for each valid task
            % and ignores bad channels.

            %% Print
            fprintf('Running function dataFilter..\n');

            %% Loop though the tasks
            for ii = 1:length(obj.validTasks)

                %% Get the task
                task = obj.validTasks{ii,1};

                %% Create filter parameters
                order = 1024;
                Wn = lfreq/(obj.(task).fs/2);
                b = fir1(order,Wn,'high');

                %% Loop through the channels and filter if they are valid
                for jj = 1:length(obj.(task).keepIdx)

                    %% Filter channel if it is valid
                    if obj.(task).keepIdx(jj,1) == 1
                        obj.(task).filteredData(jj,:) = filtfilt(b,1,obj.(task).LFP_data(jj,:)');
                    end
                end

                %% Add high pass frequency to the object
                obj.(task).lfreq = lfreq;
            end
        end

        %% Bipolar referencing
        function bipolarReference(obj)

            %% Notes
            % This function loops through good channels for each valid task
            % and performs bipolar referncing using good adjacent contacts.

            %% Print
            fprintf('Running function bipolarReference..\n');

            %% Loop through valid tasks
            for ii = 1:length(obj.validTasks)

                %% Get the task
                task = obj.validTasks{ii,1};

                %% Grab the data associated with good channels
                jj = 1;
                temp = NaN(length(obj.(task).keepIdx),length(obj.(task).filteredData));
                ch_name = cell(length(obj.(task).keepIdx),1);
                while jj <= length(obj.(task).keepIdx)
                    if obj.(task).keepIdx(jj,1) == 1
                        temp(jj,:) = obj.(task).filteredData(jj,:);
                        ch_name{jj,:} = append(obj.(task).chanSname{jj,1},num2str(jj));
                        jj = jj + 1;
                    else
                        next_idx = find(obj.(task).keepIdx(jj+1,1) == 1,1,"first");
                        if ~isempty(next_idx)
                            jj = jj + next_idx;
                            temp(jj,:) = obj.(task).filteredData(jj,:);
                            ch_name{jj,1} = append(obj.(task).chanSname{jj,1},num2str(jj));
                        end
                    end
                end

                %% Find channels that are NaN
                I = any(isnan(temp),2);
                ch_name = ch_name(~I,1);
                temp = temp(~I,:);

                %% Do bipolar referencing
                for jj = 1:length(ch_name)
                    if jj < length(ch_name)
                        obj.(task).referencedData(jj,:) =  temp(jj,:) - temp(jj+1,:);
                        obj.(task).referencedChName{jj,:} = append(ch_name{jj},'-',ch_name{jj+1,1});
                    elseif jj == length(ch_name)
                        obj.(task).referencedData(jj,:) = temp(jj,:) - temp(1,:);
                        obj.(task).referencedChName{jj,1} = append(ch_name{jj,1},'-',ch_name{1,1});
                    end
                end
            end
        end

        %% Identify TTLs
        function identifyTTLs(obj)

            %% Notes
            % This function identifies the LFP indices for cross onset
            % (11), clip on (1), clip off (2), and question (3) phases of
            % the task.  LFP indices for clip off and response are
            % calculated based on the behavior tables for each task.

            %% Print
            fprintf('Running function identifyTTLs..\n');

            %% Loop through valid tasks
            for ii = 1:length(obj.validTasks)

                %% Get the task
                task = obj.validTasks{ii,1};

                %% Manually exclude bad trials
                if strcmp(obj.subID,'MW27') & strcmp(task,'timeDiscrim')
                    obj.(task).eventIDs(128) = [];
                    obj.(task).eventTimes(128) = [];
                end

                %% Sort the TTLs of interest
                obj.(task).Indices.crossOnIdx = obj.(task).LFPIdx(obj.(task).eventIDs == 11);
                obj.(task).Indices.clipOnIdx = obj.(task).LFPIdx(obj.(task).eventIDs == 1);
                obj.(task).Indices.clipOffIdx = obj.(task).LFPIdx(obj.(task).eventIDs == 2);
                obj.(task).Indices.questionIdx = obj.(task).LFPIdx(obj.(task).eventIDs == 3);

                %% Find the resp field
                fieldNames = fieldnames(obj.(task).behave_info);
                for jj = 1:length(fieldNames)
                    if contains(fieldNames{jj,1},'respMat')
                        matchedField = fieldNames{jj,1};
                    end
                end

                %% Remove bad trial
                if strcmp(obj.subID,'MW26') & strcmp(task,'sceneRecog')
                    obj.sceneRecog.behave_info.respMat_SceneRecog(127) = []; % last trial was not completed
                    obj.(task).Indices.crossOnIdx(127) = [];
                    obj.(task).Indices.clipOnIdx(127) = [];
                    obj.(task).Indices.clipOffIdx(127) = [];
                    obj.(task).Indices.questionIdx(127) = [];
                end

                %% Create artificial indicies for crossOff and response
                cross = zeros(length(obj.(task).behave_info.(matchedField)),1);
                response = zeros(length(obj.(task).behave_info.(matchedField)),1);
                for z = 1:length((obj.(task).behave_info.(matchedField)))
                    cross(z,1) = obj.(task).behave_info.(matchedField)(z).CrossEnd - obj.(task).behave_info.(matchedField)(z).CrossStart;
                    response(z,1) = obj.(task).behave_info.(matchedField)(z).respTime - obj.(task).behave_info.(matchedField)(z).QuesStart;
                end

                %% Convert from time to samples and add to object
                obj.(task).Indices.crossOffIdx = obj.(task).Indices.crossOnIdx + (floor(obj.(task).fs .* cross));
                obj.(task).Indices.responseIdx = obj.(task).Indices.questionIdx + (floor(obj.(task).fs .* response));

            end
        end

        %% Identify boundaries
        function identifyBoundaries(obj)

            %% Notes

            %% Print
            fprintf('Running function identifyBoundaries..\n');

            %% Loop through valid tasks
            for ii = 1:length(obj.validTasks)
                
                %% Grab the task
                task = obj.validTasks{ii,1};

                %% Find the resp field
                fieldNames = fieldnames(obj.(task).behave_info);
                for jj = 1:length(fieldNames)
                    if contains(fieldNames{jj,1},'respMat')
                        matchedField = fieldNames{jj,1};
                    end
                end

                %% Identify clip field name
                x = fieldnames(obj.(task).behave_info.(matchedField));

                %% Loop through behavior info to get the task boundary type
                for jj = 1:length(obj.(task).behave_info.(matchedField))
                    obj.(task).boundary{jj,1} = extractBefore(obj.(task).behave_info.(matchedField)(jj).(x{4,1}),'_');
                end
            end
        end

        %% Allocate data
        function allocateData(obj)

            %% Notes
            % This function allocates the data based of phase of the valid
            % task (fixation, presentation, and question), as well as break
            % them up into individual epochs for SPOOOFing later on.

            %% Print
            fprintf('Running function allocateData..\n');

            %% Loop through valid tasks
            for ii = 1:length(obj.validTasks)

                %% Grab the task
                task = obj.validTasks{ii,1};

                %% Create anonymous function handle
                extractSegments = @(s,e) obj.(task).referencedData(:,s:e-1);

                %% Grab data segments
                obj.(task).Data.fixation = arrayfun(@(i) extractSegments(obj.(task).Indices.crossOnIdx(i),obj.(task).Indices.crossOffIdx(i)), 1:length(obj.(task).boundary), 'UniformOutput', false )';
                obj.(task).Data.presentation = arrayfun(@(i) extractSegments(obj.(task).Indices.clipOnIdx(i),obj.(task).Indices.clipOffIdx(i)), 1:length(obj.(task).boundary), 'UniformOutput', false )';
                obj.(task).Data.question = arrayfun(@(i) extractSegments(obj.(task).Indices.questionIdx(i),obj.(task).Indices.responseIdx(i)), 1:length(obj.(task).boundary), 'UniformOutput', false )';
                obj.(task).Data.epoch = arrayfun(@(i) extractSegments(obj.(task).Indices.crossOnIdx(i),obj.(task).Indices.responseIdx(i)), 1:length(obj.(task).boundary), 'UniformOutput', false )';

            end
        end

        %% Baseline correction
        function baselineCorrection(obj)

            %% Notes
            % This function calculates the average voltage of each channel
            % during the fixation phase and subtracts it from the voltage
            % for each channel during the presentation and question phases.

            %% Print
            fprintf('Running function baselineCorrection..\n');

            %% Loop through valid tasks
            for ii = 1:length(obj.validTasks)

                %% Grab the task
                task = obj.validTasks{ii,1};

                %% Loop through trials
                for jj = 1:length(obj.(task).Data.fixation)

                    %% Calculate channel means
                    ch_mean = mean(obj.(task).Data.fixation{jj,1},2);

                    %% Subtract from fixation
                    obj.(task).Data.fixation{jj,1} = obj.(task).Data.fixation{jj,1} - ch_mean;

                    %% Subtract from presentation
                    obj.(task).Data.presentation{jj,1} = obj.(task).Data.presentation{jj,1} - ch_mean;

                    %% Subtract from question
                    obj.(task).Data.question{jj,1} = obj.(task).Data.question{jj,1} - ch_mean;

                    %% Subtract from entire epoch
                    obj.(task).Data.epoch{jj,1} = obj.(task).Data.epoch{jj,1} - ch_mean;

                end
            end
        end

        %% Plot phases
        function plotPhases(obj)

            %% Notes
            % This function plots the phases of each task (visual check).

            %% Loop through valid tasks
            for ii = 1:length(obj.validTasks)

                %% Grab the task
                task = obj.validTasks{ii,1};

                %% Create title
                title_str = append(obj.subID,' ','Task',' ',task);

                %% Plot
                figure(ii)
                fig = gca;
                fig.XAxis.FontSize = 24;
                fig.XAxis.FontWeight = 'bold';
                fig.YAxis.FontSize = 24;
                fig.YAxis.FontWeight = 'bold';
                fig.XLim = ([min(obj.(task).LFP_converted_timestamps) max(obj.(task).LFP_converted_timestamps)]);
                fig.YLim = ([1.1*(min(min(obj.(task).referencedData))) 1.1*max(max(obj.(task).referencedData))]);
                title(title_str,'FontWeight','bold','FontSize',28);
                hold on
                for z = 1:size(obj.(task).referencedData,1)
                    plot(obj.(task).LFP_converted_timestamps,obj.(task).referencedData(z,:))
                end

                %% Figure properties
                xlabel('Time (sec)','FontWeight','bold','FontSize',24);
                ylabel('Amplitude (V)','FontWeight','bold','FontSize',24);

                %% Bands for cross presentation
                crossOnIdx = (obj.(task).Indices.crossOnIdx + obj.(task).correction_factor)./obj.(task).fs;
                crossOffIdx = (obj.(task).Indices.crossOffIdx + obj.(task).correction_factor)./obj.(task).fs;
                bands = [crossOnIdx,crossOffIdx];
                xp = [bands fliplr(bands)];
                yp = ([[1;1]*1.1*min(ylim); [1;1]*1.1*max(ylim)]*ones(1,size(bands,1))).';
                for k = 1:size(bands,1)
                    patch(xp(k,:),yp(k,:),[1 0 0],'FaceAlpha',0.1,'EdgeColor','none')
                end

                %% Bands for clip presentation
                clipOnIdx = (obj.(task).Indices.clipOnIdx + obj.(task).correction_factor)./obj.(task).fs;
                clipOffIdx = (obj.(task).Indices.clipOffIdx + obj.(task).correction_factor)./obj.(task).fs;
                bands = [clipOnIdx,clipOffIdx];
                xp = [bands fliplr(bands)];
                for k = 1:size(bands,1)
                    patch(xp(k,:),yp(k,:),[0 1 0],'FaceAlpha',0.1,'EdgeColor','none')
                end

                %% Bands for response
                questionIdx = (obj.(task).Indices.questionIdx + obj.(task).correction_factor)./obj.(task).fs;
                responseIdx = (obj.(task).Indices.responseIdx + obj.(task).correction_factor)./obj.(task).fs;
                bands = [questionIdx,responseIdx];
                xp = [bands fliplr(bands)];
                for k = 1:size(bands,1)
                    patch(xp(k,:),yp(k,:),[0 0 1],'FaceAlpha',0.1,'EdgeColor','none')
                end

            end
        end

        %% Compute Power
        function computePower(obj,freq)

            %% Notes
            % This function computes the power of all phases for each task
            % in the object.  In the event there is not enough data to
            % capture the frequency range of interest, NaNs are added to
            % the appropriate field of the object.

            %% Print
            fprintf('Running function computePower..\n');

            %% Set valid phases
            obj.validPhases = {"fixation";"presentation";"question"};

            %% Loop through valid tasks
            for ii = 1:length(obj.validTasks)

                %% Get the task
                task = obj.validTasks{ii,1};

                %% The the frequencies
                obj.(task).powerFreq = freq;

                %% Loop through valid phases
                for jj = 1:length(obj.validPhases)

                    %% Loop through trials
                    for kk = 1:length(obj.(task).Data.(obj.validPhases{jj,1}))

                        %% Grab the trial data
                        tmp = obj.(task).Data.(obj.validPhases{jj,1}){kk,1};

                        %% Loop through channels and compute CWT
                        for ll = 1:size(tmp,1)

                            %% Create field in the object
                            powerField = append(obj.validPhases{jj,1},'Power');

                            %% Ignore NaN
                            if ~isnan(tmp(ll,:)) & length(tmp) > 200

                                %% Compute power
                                [wt,f] = cwt(tmp(ll,:),obj.(task).fs);

                                %% Grab power for frequency range of interest
                                I = f >= freq(1,1) & f <= freq(2,1);
                                obj.(task).Power.(powerField){kk,1}(ll,:) = mean(abs(wt(I,:)),1);

                            else
                                obj.(task).Power.(powerField){kk,1}(ll,:) = NaN;
                            end

                        end
                    end
                end
            end
        end

        %% Compute statistics for power differences
        function computeStats(obj)

            %% Notes
            % This function computes the statistics for power changes
            % between the presentation and question phases for each valid
            % task.  Stats are corrected for multiple comparisons using
            % Benjamini-Hochberg, and the statistically significant
            % channels and their associated p-values are added to the
            % respective fields of the object.

            %% Print
            fprintf('Running function computeStats..\n');

            %% Loop over task
            for ii = 1:length(obj.validTasks)

                %% Grab the task
                task = obj.validTasks{ii,1};

                %% Loop over phases
                for jj = 1:length(obj.validPhases)

                    %% Grab the phase
                    phase = obj.validPhases{jj,1};
                    phasePower = append(phase,'Power');

                    %% Loop through trials to calculate average power
                    for kk = 1:length(obj.(task).Power.(phasePower))

                        %% Create field name
                        ave_power = append('ave_',phasePower);

                        %% Grab the data
                        tmp = obj.(task).Power.(phasePower){kk,1};

                        %% Check if any data isnan
                        if ~any(isnan(tmp))

                            %% Loop through to calculate average power
                            for ll = 1:size(tmp,1)
                                obj.(task).Power.(ave_power)(ll,kk) = mean(obj.(task).Power.(phasePower){kk,1}(ll,:));
                            end
                        end
                    end
                end

                %% Compare power in presentation vs question phase
                for jj = 1:size(obj.(task).Power.ave_presentationPower,1)
                    obj.(task).Power.powerStats(jj,1) = ranksum(obj.(task).Power.ave_presentationPower(jj,:), obj.(task).Power.ave_questionPower(jj,:));
                end

                %% Correct for false discovery using Benjamini-Hochberg
                [val,idx] = sort(obj.(task).Power.powerStats,'ascend');
                adjusted_p = zeros(length(val),1);
                corrected = NaN(length(val),1);
                for z = 1:length(val)
                    adjusted_p(z,1) = 0.05/(length(val)+1-z);
                    if val(z,1) <= adjusted_p(z,1)
                        corrected(z,1) = 1;
                    else
                        corrected(z,1) = 0;
                    end
                end

                %% Identify significant channels
                corrected = logical(corrected);
                new_names_corrected = obj.(task).referencedChName(idx);
                obj.(task).Power.powerProb = val(corrected);
                obj.(task).Power.powerSigChannels = new_names_corrected(corrected);

            end
        end

        %% Compute functional connectivity
        function functionalConnectivity(obj,method)

            %% Notes
            % Granger FC uses lowest value for power frequency used earlier
            % to calculate number of lags.
            %
            % This function calculates the functional connectivity of each
            % valid phase for each task 

            %% Print
            fprintf('Running function functionalConnectivity..\n');

            %% Loop through valid tasks
            for ii = 1:length(obj.validTasks)

                %% Grab the task
                task = obj.validTasks{ii,1};

                %% Loop through phases
                for jj = 1:length(obj.validPhases)

                    %% Grab the phase
                    phase = obj.validPhases{jj,1};

                    %% Add method to object
                    obj.(task).Connectivity.fcMethod = method;

                    %% Loop through trials
                    for kk = 1:length(obj.(task).Data.(phase))

                        %% Grab the data
                        tmp = obj.(task).Data.(phase){kk,1};

                        %% Compute FC based on method
                        if strcmp(method,'granger')

                            %% Calculate number of samples
                            lags = ceil(1/min(obj.(task).powerFreq) * obj.(task).fs);

                            %% Preallocate
                            fc = append(phase,'FC');
                            obj.(task).Connectivity.(fc){kk,1} = zeros(size(tmp,1),size(tmp,1));

                            %% Skip trial if not enough data
                            if length(tmp) < 4*lags
                                obj.(task).Connectivity.(fc){kk,1} = NaN(size(tmp,1),size(tmp,1)); % or zeros
                                continue
                            end

                            %% Add lags to the object
                            obj.(task).Connectivity.numLags = lags;

                            %% Loop over channels
                            for ll = 1:size(tmp,1)
                                for mm = 1:size(tmp,1)
                                    if ll ~= mm

                                        %% Grab the data
                                        sig_reduced = tmp(ll,:);
                                        sig_full = [tmp(ll,:);tmp(mm,:)];

                                        %% Calculate residuals
                                        [R_reduced,~] = AutoregressiveProcess(sig_reduced,lags);
                                        [R_full,~] = AutoregressiveProcess(sig_full,lags);

                                        %% Calculate variances
                                        var_reduced = var(R_reduced(1,:),0,2);
                                        var_full = var(R_full(1,:),0,2);

                                        %% Write FC values
                                        obj.(task).Connectivity.(fc){kk,1}(ll,mm) = log(var_reduced/var_full);

                                        %% Clear variables
                                        clear sig_reduced sig_full R_reduced R_full var_reduced var_full;

                                    end
                                end
                            end
                        elseif strcmp(method,'coherence')

                            %% Preallocate
                            fc = append(phase,'FC');
                            obj.(task).Connectivity.(fc){kk,1} = zeros(size(tmp,1),size(tmp,1));

                            %% Loop over channels
                            for ll = 1:size(tmp,1)
                                for mm = 1:size(tmp,1)
                                    if ll ~= mm

                                        %% Calculate magnitude squared coherence
                                        [x,f] = mscohere(tmp(ll,:),tmp(mm,:),[],[],[],obj.(task).fs);

                                        %% Grab the values in the frequency range
                                        J = f >= min(obj.(task).powerFreq) & f <= max(obj.(task).powerFreq);

                                        %% Write the mean values
                                        obj.(task).Connectivity.(fc){kk,1}(ll,mm) = mean(x(J));
                                    elseif ll == mm
                                        fc = append(phase,'FC');
                                        obj.(task).(fc){kk,1}(ll,mm) = 0;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end

        %% Analyze FC
        function analyzeFC(obj,boundary)

            %% Notes
            % This function calculates the statistical significance of the
            % difference between functiional connectivity in the
            % presentation phase and the question phase for all valid
            % tasks.  Statistics and significant channel combinations are
            % added to the object.

            %% Print
            fprintf('Running function analyzeFC..\n');

            %% Loop through the tasks
            for ii = 1:length(obj.validTasks)

                %% Get the task
                task = obj.validTasks{ii,1};

                %% Compare FC between presentation and question phases
                if isfield(obj.(task).Connectivity,'presentationFC') && isfield(obj.(task).Connectivity,'questionFC')

                    %% Identify boundary of interest
                    I = strcmp(obj.(task).boundary,boundary);
                    presentation = obj.(task).Connectivity.presentationFC(I);
                    question = obj.(task).Connectivity.questionFC(I);

                    %% Calculate mean FC for the phases
                    pres_fc = zeros(length(obj.(task).referencedChName),length(obj.(task).referencedChName));
                    pres_count = 0;
                    for jj = 1:length(presentation)
                        if ~any(isnan(presentation{jj,1}(:,:)))
                            pres_count = pres_count + 1;
                            pres_fc = pres_fc + presentation{jj,1}(:,:);
                        else
                            continue
                        end
                    end
                    pres_fc = pres_fc ./ pres_count;

                    ques_fc = zeros(length(obj.(task).referencedChName),length(obj.(task).referencedChName));
                    ques_count = 0;
                    for jj = 1:length(question)
                        if ~any(isnan(question{jj,1}(:,:)))
                            ques_count = ques_count + 1;
                            ques_fc = ques_fc + question{jj,1}(:,:);
                        else
                            continue
                        end
                    end
                    ques_fc = ques_fc ./ ques_count;

                    fc_presentation_max = max(max(pres_fc));
                    fc_question_max = max(max(ques_fc));
                    y = max(fc_question_max,fc_presentation_max);

                    fc_presentation_min = min(min(pres_fc));
                    fc_question_min = min(min(ques_fc));
                    x = min(fc_presentation_min,fc_question_min);

                    title1_string = append(obj.subID,' ',task,' Presentation FC');
                    title2_string = append(obj.subID,' ',task,' Question FC');

                    figure
                    imagesc(pres_fc);
                    title(title1_string,'FontSize',24,'FontWeight','bold');
                    fig = gca;
                    fig.XAxis.FontSize = 24;
                    fig.XAxis.FontWeight = 'bold';
                    xlabel('Channel (From)','FontWeight','bold','FontSize',24);
                    fig.YAxis.FontSize = 24;
                    fig.YAxis.FontWeight = 'bold';
                    ylabel('Channel (To)','FontWeight','bold','FontSize',24);
                    clim manual;
                    clim([x y]);
                    cb = colorbar;
                    cb.Label.String = 'F-Value';
                    cb.FontSize = 24;
                    cb.FontWeight = 'bold';cb.FontWeight = 'bold';

                    figure
                    imagesc(ques_fc);
                    title(title2_string,'FontSize',24,'FontWeight','bold');
                    fig = gca;
                    fig.XAxis.FontSize = 24;
                    fig.XAxis.FontWeight = 'bold';
                    xlabel('Channel (From)','FontWeight','bold','FontSize',24);
                    fig.YAxis.FontSize = 24;
                    fig.YAxis.FontWeight = 'bold';
                    ylabel('Channel (To)','FontWeight','bold','FontSize',24);
                    clim manual;
                    clim([x y]);
                    cb = colorbar;
                    cb.Label.String = 'F-Value';
                    cb.FontSize = 24;
                    cb.FontWeight = 'bold';cb.FontWeight = 'bold';

                    %% Calculate probability values for FC
                    p = ones(length(obj.(task).referencedChName)*length(obj.(task).referencedChName),1);
                    ch_names = cell(length(obj.(task).referencedChName)*length(obj.(task).referencedChName),1);
                    ch_ref = cell(length(obj.(task).referencedChName),length(obj.(task).referencedChName));
                    count = 1;
                    for jj = 1:length(obj.(task).referencedChName)
                        for kk = 1:length(obj.(task).referencedChName)
                            pres_val = NaN(length(presentation),1);
                            ques_val = NaN(length(question),1);
                            ch_ref{jj,kk} = append(obj.(task).referencedChName{jj,1},'\',obj.(task).referencedChName{kk,1});
                            if jj ~= kk
                                for ll = 1:length(question)
                                    if ~any(isnan(question{ll,1}))
                                        pres_val(ll,1) = presentation{ll,1}(jj,kk);
                                        ques_val(ll,1) = question{ll,1}(jj,kk);
                                    else
                                        continue
                                    end
                                end
                                [~,p(count,1)] = ttest(pres_val,ques_val);
                                ch_names{count,1} = append(obj.(task).referencedChName{jj,1},'\',obj.(task).referencedChName{kk,1});
                                count = count + 1;
                            else
                                ch_names{count,1} = append(obj.(task).referencedChName{jj,1},'\',obj.(task).referencedChName{kk,1});
                                count = count + 1;
                                continue
                            end
                        end
                    end

                    %% Correct for false discovery using Benjamini-Hochberg
                    [val,idx] = sort(p,'ascend');
                    adjusted_p = zeros(length(val),1);
                    corrected = NaN(length(val),1);
                    for z = 1:length(val)
                        adjusted_p(z,1) = 0.05/(length(val)+1-z);
                        if val(z,1) <= adjusted_p(z,1)
                            corrected(z,1) = 1;
                        else
                            corrected(z,1) = 0;
                        end
                    end

                    %% Identify significant channels
                    corrected = logical(corrected);
                    new_names_corrected = ch_names(idx);
                    obj.(task).Connectivity.fcprob = val(corrected);
                    obj.(task).Connectivity.ch_names_corrected = new_names_corrected(corrected);

                end
            end
        end

        %% SPOOOF
        function SPOOOF(obj,l_freq,h_freq)

            %% Notes
            % This function performs SPRiNT-FOOOF (SPOOOFing) on each
            % channel for each epoch for each task, and adds the result to
            % the object.

            %% Print
            fprintf('Running function SPOOOF..\n');

            %% Loop through tasks
            for ii = 1:length(obj.validTasks)

                %% Grab the task
                task = obj.validTasks{ii,1};

                %% Preallocate
                opt = struct();

                %% Inputs
                % STFT options
                opt.sfreq = obj.(task).fs; % input sampling rate
                opt.WinLength = 0.5; % STFT window length; originally 1 (one second) (Lisa had 0.5)
                opt.WinOverlap = 50; % Overlap between sliding windows (in %)
                opt.WinAverage = 2; % Number of sliding windows averaged by time point

                % specparams options
                opt.freq_range = [l_freq h_freq];
                opt.peak_width_limits   = [1.5 7];
                opt.max_peaks           = 8;
                opt.min_peak_height     = 6 / 10; % convert from dB to B
                opt.aperiodic_mode      = 'knee'; % alternative: fixed
                opt.peak_threshold      = 2; % 2 std dev: parameter for interface simplification, Lisa had 1.5

                % MATLAB-only options
                opt.peak_type           = 'gaussian'; % alternative: cauchy
                opt.proximity_threshold = 2;
                opt.guess_weight        = 'none';
                opt.thresh_after        = true;   % Threshold after fitting, always selected for Matlab

                % (mirrors the Python FOOOF closest by removing peaks
                % that do not satisfy a user's predetermined conditions)
                % only used in the absence of the

                if license('test','optimization_toolbox') % check for optimization toolbox
                    opt.hOT = 1;
                    disp('Using constrained optimization, Guess Weight ignored.')
                else
                    opt.hOT = 0;
                    disp('Using unconstrained optimization, with Guess Weights.')
                end

                opt.rmoutliers          = 'yes';
                opt.maxfreq             = 2.5;
                opt.maxtime             = 6;
                opt.minnear             = 3;

                Freqs = 0:1/opt.WinLength:opt.sfreq/2;

                %% Loop over epochs
                for jj = 1:length(obj.(task).Data.epoch)

                    %% Print epoch
                    fprintf('Analyzing epoch %s..\n',num2str(jj))

                    %% Loop over channe;s
                    for kk = 1:length(obj.(task).referencedChName)

                        %% Get temp data
                        tempChan = obj.(task).Data.epoch{jj,1}(kk,:);

                        %% Create structure
                        channel = struct('data',[],'peaks',[],'aperiodics',[],'stats',[]);

                        %% Compute short-time Fourier transform
                        [TF, ts] = SPRiNT_stft(tempChan,opt);
                        outputStruct = struct('opts',opt,'freqs',Freqs,'channel',channel);

                        %% Parameterize STFTs
                        s_data = SPRiNT_specparam_matlab(TF,outputStruct.freqs,outputStruct.opts,ts);
                        s_data.SPRiNT.SPRiNT_models = squeeze(s_data.SPRiNT.SPRiNT_models);
                        s_data.SPRiNT.peak_models = squeeze(s_data.SPRiNT.peak_models);
                        s_data.SPRiNT.aperiodic_models = squeeze(s_data.SPRiNT.aperiodic_models);

                        sprintOUT = s_data.SPRiNT;
                        obj.(task).SPRiNTout{jj,1}{kk,1} = sprintOUT;

                    end
                end
            end
        end

        %% Next function

    end
end

%% SPRiNT_stft
function [TF, ts] = SPRiNT_stft(F,opts)
    %% Notes
    % SPRiNT_stft: Compute a locally averaged short-time Fourier transform (for
    % use in SPRiNT)
    %
    % Segments of this function were adapted from the Brainstorm software package:
    % https://neuroimage.usc.edu/brainstorm
    % Tadel et al. (2011)
    %
    % Copyright (c)2000-2020 University of Southern California & McGill University
    % This software is distributed under the terms of the GNU General Public License
    % as published by the Free Software Foundation. Further details on the GPLv3
    % license can be found at http://www.gnu.org/copyleft/gpl.html.
    %
    % Author: Luc Wilson (2022)

    %% Set parameters
    sfreq = opts.sfreq; % sample rate, in Hertz
    WinLength = opts.WinLength; % window length, in seconds
    WinOverlap = opts.WinOverlap; % window overlap, in percent
    avgWin = opts.WinAverage; % number of windows being averaged per PSD
    nTime = size(F,2); % number of time points

    %% Windowing
    Lwin  = round(WinLength * sfreq); % number of data points in windows
    Loverlap = round(Lwin * WinOverlap / 100); % number of data points in overlap

    %% Check if window is too small
    if (Lwin < 50)
        return;
        % Window is larger than the data
    elseif (Lwin > nTime)
        Lwin = size(F,2);
        Lwin = Lwin - mod(Lwin,2); % Make sure the number of samples is even
        Loverlap = 0;
        Nwin = 1;
    else
        Lwin = Lwin - mod(Lwin,2); % Make sure the number of samples is even
        Nwin = floor((nTime - Loverlap) ./ (Lwin - Loverlap));
    end

    %% Other parameters
    % Next power of two from length of signal
    NFFT = Lwin; % No zero-padding: Nfft = Ntime
    
    % Positive frequency bins spanned by FFT
    FreqVector = sfreq / 2 * linspace(0,1,NFFT/2+1);

    % Determine Hann window shape/power
    Win = hann(Lwin)';
    WinNoisePowerGain = sum(Win.^2);

    %% Initialize STFT, time matrices
    ts = nan(Nwin-(avgWin-1),1);
    TF = nan(size(F,1), Nwin-(avgWin-1), size(FreqVector,2));
    TFtmp = nan(size(F,1), avgWin, size(FreqVector,2));

    %% Calculate FFT for each window
    TFfull = zeros(size(F,1),Nwin,size(FreqVector,2));
    for iWin = 1:Nwin

        % Build indices
        iTimes = (1:Lwin) + (iWin-1)*(Lwin - Loverlap);
        center_time = floor(median((iTimes-(avgWin-1)./2*(Lwin - Loverlap))))./sfreq;

        % Select indices
        Fwin = F(:,iTimes);

        % No need to enforce removing DC component (0 frequency).
        Fwin = Fwin - mean(Fwin,2);

        % Apply Hann window to signal
        Fwin = Fwin .* Win;

        % Compute FFT
        Ffft = fft(Fwin, NFFT, 2);

        % One-sided spectrum (keep only first half)
        % (x2 to recover full power from negative frequencies)
        TFwin = Ffft(:,1:NFFT/2+1) * sqrt(2 ./ (sfreq * WinNoisePowerGain));

        % x2 doesn't apply to DC and Nyquist
        TFwin(:, [1,end]) = TFwin(:, [1,end]) ./ sqrt(2);

        % Permute dimensions: time and frequency
        TFwin = permute(TFwin, [1 3 2]);

        % Convert to power
        TFwin = abs(TFwin).^2;
        TFfull(:,iWin,:) = TFwin;
        TFtmp(:,mod(iWin,avgWin)+1,:) = TFwin;
        if isnan(TFtmp(1,1,1))
            continue % Do not record anything until transient is gone
        else
            % Save STFTs for window
            TF(:,iWin-(avgWin-1),:) = mean(TFtmp,2);
            ts(iWin-(avgWin-1)) = center_time;
        end
    end
end

%% Step 2: parameterize spectrograms spectra
function [s_data] = SPRiNT_specparam_matlab(TF, fs, opt, ts)

    %% Notes
    % SPRiNT_specparam_matlab: Compute time-resolved specparam models for
    % short-time Fourier transformed signals.
    %
    % The spectral parameterization algorithm used herein (specparam) can be
    % cited as:
    %   Donoghue, T., Haller, M., Peterson, E.J. et al., Parameterizing neural
    %   power spectra into periodic and aperiodic components. Nat Neurosci 23,
    %   1655–1665 (2020). https://doi.org/10.1038/s41593-020-00744-x
    %
    % Author: Luc Wilson (2022)

    %% Parameters
    fMask = logical(round(fs.*10)./10 >= round(opt.freq_range(1).*10)./10 & (round(fs.*10)./10 <= round(opt.freq_range(2).*10)./10));
    fs = fs(fMask);
    s_data.Freqs = fs;
    nChan = size(TF,1);
    nTimes = size(TF,2);

    %% Adjust TF plots to only include modelled frequencies
    TF = TF(:,:,fMask);

    %% Initialize FOOOF structs
    channel(nChan) = struct();
    SPRiNT = struct('options',opt,'freqs',fs,'channel',channel,'SPRiNT_models',nan(size(TF)),'peak_models',nan(size(TF)),'aperiodic_models',nan(size(TF)));

    %% Iterate across channels
    for chan = 1:nChan
        channel(chan).data(nTimes) = struct(...
            'time',             [],...
            'aperiodic_params', [],...
            'peak_params',      [],...
            'peak_types',       '',...
            'ap_fit',           [],...
            'fooofed_spectrum', [],...
            'power_spectrum',   [],...
            'peak_fit',         [],...
            'error',            [],...
            'r_squared',        []);
        channel(chan).peaks(nTimes*opt.max_peaks) = struct(...
            'time',             [],...
            'center_frequency', [],...
            'amplitude',        [],...
            'st_dev',           []);
        channel(chan).aperiodics(nTimes) = struct(...
            'time',             [],...
            'offset',           [],...
            'exponent',         []);
        channel(chan).stats(nTimes) = struct(...
            'MSE',              [],...
            'r_squared',        [],...
            'frequency_wise_error', []);
        spec = log10(squeeze(TF(chan,:,:))); % extract log spectra for a given channel

        % Iterate across time
        i = 1; % For peak extraction
        ag = -(spec(1,end)-spec(1,1))./log10(fs(end)./fs(1)); % aperiodic guess initialization
        for time = 1:nTimes

            %% Fit aperiodic
            aperiodic_pars = robust_ap_fit(fs, spec(time,:), opt.aperiodic_mode, ag);

            %% Remove aperiodic
            flat_spec = flatten_spectrum(fs, spec(time,:), aperiodic_pars, opt.aperiodic_mode);

            %% Fit peaks
            [peak_pars, peak_function] = fit_peaks(fs, flat_spec, opt.max_peaks, opt.peak_threshold, opt.min_peak_height, ...
                opt.peak_width_limits/2, opt.proximity_threshold, opt.peak_type, opt.guess_weight,opt.hOT);

            %% Check thresholding requirements are met for unbounded optimization
            if opt.thresh_after && ~opt.hOT
                peak_pars(peak_pars(:,2) < opt.min_peak_height,:)     = []; % remove peaks shorter than limit
                peak_pars(peak_pars(:,3) < opt.peak_width_limits(1)/2,:)  = []; % remove peaks narrower than limit
                peak_pars(peak_pars(:,3) > opt.peak_width_limits(2)/2,:)  = []; % remove peaks broader than limit
                peak_pars = drop_peak_cf(peak_pars, opt.proximity_threshold, opt.freq_range); % remove peaks outside frequency limits
                peak_pars(peak_pars(:,1) < 0,:) = []; % remove peaks with a centre frequency less than zero (bypass drop_peak_cf)
                peak_pars = drop_peak_overlap(peak_pars, opt.proximity_threshold); % remove smallest of two peaks fit too closely
            end

            %% Refit aperiodic
            aperiodic = spec(time,:);
            for peak = 1:size(peak_pars,1)
                aperiodic = aperiodic - peak_function(fs,peak_pars(peak,1), peak_pars(peak,2), peak_pars(peak,3));
            end
            aperiodic_pars = simple_ap_fit(fs, aperiodic, opt.aperiodic_mode, aperiodic_pars(end));
            ag = aperiodic_pars(end); % save exponent estimate for next iteration

            %% Generate model fit
            ap_fit = gen_aperiodic(fs, aperiodic_pars, opt.aperiodic_mode);
            model_fit = ap_fit;
            for peak = 1:size(peak_pars,1)
                model_fit = model_fit + peak_function(fs,peak_pars(peak,1),...
                    peak_pars(peak,2),peak_pars(peak,3));
            end

            %% Calculate model error
            MSE = sum((spec(time,:) - model_fit).^2)/length(model_fit);
            rsq_tmp = corrcoef(spec(time,:),model_fit).^2;

            %% Return FOOOF results
            aperiodic_pars(2) = abs(aperiodic_pars(2));
            channel(chan).data(time).time                = ts(time);
            channel(chan).data(time).aperiodic_params    = aperiodic_pars;
            channel(chan).data(time).peak_params         = peak_pars;
            channel(chan).data(time).peak_types          = func2str(peak_function);
            channel(chan).data(time).ap_fit              = 10.^ap_fit;
            aperiodic_models(chan,time,:)                = 10.^ap_fit;
            channel(chan).data(time).fooofed_spectrum    = 10.^model_fit;
            SPRiNT_models(chan,time,:)                   = 10.^model_fit;
            channel(chan).data(time).power_spectrum   	 = 10.^spec(time,:);
            channel(chan).data(time).peak_fit            = 10.^(model_fit-ap_fit);
            peak_models(chan,time,:)                     = 10.^(model_fit-ap_fit);
            channel(chan).data(time).error               = MSE;
            channel(chan).data(time).r_squared           = rsq_tmp(2);

            %% Extract peaks
            if ~isempty(peak_pars) & any(peak_pars)
                for p = 1:size(peak_pars,1)
                    channel(chan).peaks(i).time = ts(time);
                    channel(chan).peaks(i).center_frequency = peak_pars(p,1);
                    channel(chan).peaks(i).amplitude = peak_pars(p,2);
                    channel(chan).peaks(i).st_dev = peak_pars(p,3);
                    i = i + 1;
                end
            end

            %% Extract aperiodic
            channel(chan).aperiodics(time).time = ts(time);
            channel(chan).aperiodics(time).offset = aperiodic_pars(1);
            if length(aperiodic_pars)>2 % Legacy specparam alters order of parameters
                channel(chan).aperiodics(time).exponent = aperiodic_pars(3);
                channel(chan).aperiodics(time).knee_frequency = aperiodic_pars(2);
            else
                channel(chan).aperiodics(time).exponent = aperiodic_pars(2);
            end
            channel(chan).stats(time).MSE = MSE;
            channel(chan).stats(time).r_squared = rsq_tmp(2);
            channel(chan).stats(time).frequency_wise_error = abs(spec(time,:)-model_fit);

        end
        channel(chan).peaks(i:end) = [];
    end
    SPRiNT.channel = channel;
    SPRiNT.aperiodic_models = aperiodic_models;
    SPRiNT.SPRiNT_models = SPRiNT_models;
    SPRiNT.peak_models = peak_models;

    if strcmp(opt.rmoutliers,'yes')
        SPRiNT = remove_outliers(SPRiNT,peak_function,opt);
    end

    SPRiNT = cluster_peaks_dynamic2(SPRiNT); % Cluster peaks
    s_data.SPRiNT = SPRiNT;
end

%% Fit aperiodic
function aperiodic_params = robust_ap_fit(freqs, power_spectrum, aperiodic_mode, aperiodic_guess)

    %% Notes
    %       Fit the aperiodic component of the power spectrum robustly, ignoring outliers.
    %
    %       Parameters
    %       ----------
    %       freqs : 1xn array
    %           Frequency values for the power spectrum, in linear scale.
    %       power_spectrum : 1xn array
    %           Power values, in log10 scale.
    %       aperiodic_mode : {'fixed','knee'}
    %           Defines absence or presence of knee in aperiodic component.
    %       aperiodic_guess: double
    %           SPRiNT specific - feeds previous timepoint aperiodic slope as
    %           guess
    %
    %       Returns
    %       -------
    %       aperiodic_params : 1xn array
    %           Parameter estimates for aperiodic fit.

    %% Do a quick, initial aperiodic fit
    popt = simple_ap_fit(freqs, power_spectrum, aperiodic_mode, aperiodic_guess);
    initial_fit = gen_aperiodic(freqs, popt, aperiodic_mode);

    %% Flatten power spectrum based on initial aperiodic fit
    flatspec = power_spectrum - initial_fit;

    %% Flatten outliers (any points that drop below 0)
    flatspec(flatspec(:) < 0) = 0;

    %% Use percential threshold, in terms of # of points, to extract and re-fit
    perc_thresh = prctile(flatspec, 0.025);
    perc_mask = flatspec <= perc_thresh;
    freqs_ignore = freqs(perc_mask);
    spectrum_ignore = power_spectrum(perc_mask);

    %% Second aperiodic fit - using results of first fit as guess parameters
    options = optimset('Display', 'off', 'TolX', 1e-4, 'TolFun', 1e-6, ...
        'MaxFunEvals', 5000, 'MaxIter', 5000);
    guess_vec = popt;

    switch (aperiodic_mode)
        case 'fixed'  % no knee
            aperiodic_params = fminsearch(@error_expo_nk_function, guess_vec, options, freqs_ignore, spectrum_ignore);
        case 'knee'
            aperiodic_params = fminsearch(@error_expo_function, guess_vec, options, freqs_ignore, spectrum_ignore);
    end
end

%% Fitting algorithm
function aperiodic_params = simple_ap_fit(freqs, power_spectrum, aperiodic_mode, aperiodic_guess)

    %% Notes
    %       Fit the aperiodic component of the power spectrum.
    %
    %       Parameters
    %       ----------
    %       freqs : 1xn array
    %           Frequency values for the power spectrum, in linear scale.
    %       power_spectrum : 1xn array
    %           Power values, in log10 scale.
    %       aperiodic_mode : {'fixed','knee'}
    %           Defines absence or presence of knee in aperiodic component.
    %       aperiodic_guess: double
    %           SPRiNT specific - feeds previous timepoint aperiodic slope as
    %           guess
    %
    %       Returns
    %       -------
    %       aperiodic_params : 1xn array
    %           Parameter estimates for aperiodic fit.

    %       Set guess params for lorentzian aperiodic fit, guess params set at init

    %%
    options = optimset('Display', 'off', 'TolX', 1e-4, 'TolFun', 1e-6, ...
        'MaxFunEvals', 5000, 'MaxIter', 5000);

    switch (aperiodic_mode)
        case 'fixed'  % no knee
            guess_vec = [power_spectrum(1), aperiodic_guess];
            aperiodic_params = fminsearch(@error_expo_nk_function, guess_vec, options, freqs, power_spectrum);
        case 'knee'
            guess_vec = [power_spectrum(1),0, aperiodic_guess];
            aperiodic_params = fminsearch(@error_expo_function, guess_vec, options, freqs, power_spectrum);
    end
end

%% Flatten spectrum
function spectrum_flat = flatten_spectrum(freqs, power_spectrum, robust_aperiodic_params, aperiodic_mode)

    %% Notes
    %       Flatten the power spectrum by removing the aperiodic component.
    %
    %       Parameters
    %       ----------
    %       freqs : 1xn array
    %           Frequency values for the power spectrum, in linear scale.
    %       power_spectrum : 1xn array
    %           Power values, in log10 scale.
    %       robust_aperiodic_params : 1x2 or 1x3 array (see aperiodic_mode)
    %           Parameter estimates for aperiodic fit.
    %       aperiodic_mode : 1 or 2
    %           Defines absence or presence of knee in aperiodic component.
    %
    %       Returns
    %       -------
    %       spectrum_flat : 1xn array
    %           Flattened (aperiodic removed) power spectrum.

    %%
    spectrum_flat = power_spectrum - gen_aperiodic(freqs,robust_aperiodic_params,aperiodic_mode);

end

%% Generate aperiodic
function ap_vals = gen_aperiodic(freqs,aperiodic_params,aperiodic_mode)

    %% Notes
    %       Generate aperiodic values, from parameter definition.
    %
    %       Parameters
    %       ----------
    %       freqs : 1xn array
    %       	Frequency vector to create aperiodic component for.
    %       aperiodic_params : 1x3 array
    %           Parameters that define the aperiodic component.
    %       aperiodic_mode : {'fixed', 'knee'}
    %           Defines absence or presence of knee in aperiodic component.
    %
    %       Returns
    %       -------
    %       ap_vals : 1d array
    %           Generated aperiodic values.

    %%
    switch aperiodic_mode
        case 'fixed'  % no knee
            ap_vals = expo_nk_function(freqs,aperiodic_params);
        case 'knee'
            ap_vals = expo_function(freqs,aperiodic_params);
    end
end

%% Fit peaks
function [model_params,peak_function] = fit_peaks(freqs, flat_iter, max_n_peaks, peak_threshold, min_peak_height, gauss_std_limits, proxThresh, peakType, guess_weight, hOT)

    %% Notes
    %       Iteratively fit peaks to flattened spectrum.
    %
    %       Parameters
    %       ----------
    %       freqs : 1xn array
    %           Frequency values for the power spectrum, in linear scale.
    %       flat_iter : 1xn array
    %           Flattened (aperiodic removed) power spectrum.
    %       max_n_peaks : double
    %           Maximum number of gaussians to fit within the spectrum.
    %       peak_threshold : double
    %           Threshold (in standard deviations of noise floor) to detect a peak.
    %       min_peak_height : double
    %           Minimum height of a peak (in log10).
    %       gauss_std_limits : 1x2 double
    %           Limits to gaussian (cauchy) standard deviation (gamma) when detecting a peak.
    %       proxThresh : double
    %           Minimum distance between two peaks, in st. dev. (gamma) of peaks.
    %       peakType : {'gaussian', 'cauchy'}
    %           Which types of peaks are being fitted
    %       guess_weight : {'none', 'weak', 'strong'}
    %           Parameter to weigh initial estimates during optimization (None, Weak, or Strong)
    %       hOT : 0 or 1
    %           Defines whether to use constrained optimization, fmincon, or
    %           basic simplex, fminsearch.
    %
    %       Returns
    %       -------
    %       gaussian_params : mx3 array, where m = No. of peaks.
    %           Parameters that define the peak fit(s). Each row is a peak, as [mean, height, st. dev. (gamma)].

    %% Peak type
    switch peakType
        case 'gaussian' % gaussian only
            peak_function = @gaussian; % Identify peaks as gaussian

            % Initialize matrix of guess parameters for gaussian fitting.
            guess_params = zeros(max_n_peaks, 3);

            % Save intact flat_spectrum
            flat_spec = flat_iter;

            % Find peak: Loop through, finding a candidate peak, and fitting with a guess gaussian.
            % Stopping procedure based on either the limit on # of peaks,
            % or the relative or absolute height thresholds.
            for guess = 1:max_n_peaks

                % Find candidate peak - the maximum point of the flattened spectrum.
                max_ind = find(flat_iter == max(flat_iter));
                max_height = flat_iter(max_ind);

                % Stop searching for peaks once max_height drops below height threshold.
                if max_height <= peak_threshold * std(flat_iter)
                    break
                end

                % Set the guess parameters for gaussian fitting - mean and height.
                guess_freq = freqs(max_ind);
                guess_height = max_height;

                % Halt fitting process if candidate peak drops below minimum height.
                if guess_height <= min_peak_height
                    break
                end

                % Data-driven first guess at standard deviation
                % Find half height index on each side of the center frequency.
                half_height = 0.5 * max_height;

                le_ind = sum(flat_iter(1:max_ind) <= half_height);
                ri_ind = length(flat_iter) - sum(flat_iter(max_ind:end) <= half_height)+1;

                % Keep bandwidth estimation from the shortest side.
                % We grab shortest to avoid estimating very large std from overalapping peaks.
                % Grab the shortest side, ignoring a side if the half max was not found.
                % Note: will fail if both le & ri ind's end up as None (probably shouldn't happen).
                short_side = min(abs([le_ind,ri_ind]-max_ind));

                % Estimate std from FWHM. Calculate FWHM, converting to Hz, get guess std from FWHM
                fwhm = short_side * 2 * (freqs(2)-freqs(1));
                guess_std = fwhm / (2 * sqrt(2 * log(2)));

                % Check that guess std isn't outside preset std limits; restrict if so.
                % Note: without this, curve_fitting fails if given guess > or < bounds.
                if guess_std < gauss_std_limits(1)
                    guess_std = gauss_std_limits(1);
                end
                if guess_std > gauss_std_limits(2)
                    guess_std = gauss_std_limits(2);
                end

                % Collect guess parameters.
                guess_params(guess,:) = [guess_freq, guess_height, guess_std];

                % Subtract best-guess gaussian.
                peak_gauss = gaussian(freqs, guess_freq, guess_height, guess_std);
                flat_iter = flat_iter - peak_gauss;

            end

            % Remove unused guesses
            guess_params(guess_params(:,1) == 0,:) = [];

            % Check peaks based on edges, and on overlap
            % Drop any that violate requirements.
            guess_params = drop_peak_cf(guess_params, proxThresh, [min(freqs) max(freqs)]);
            guess_params = drop_peak_overlap(guess_params, proxThresh);

            % If there are peak guesses, fit the peaks, and sort results.
            if ~isempty(guess_params)
                model_params = fit_peak_guess(guess_params, freqs, flat_spec, 1, guess_weight, gauss_std_limits,hOT);
            else
                model_params = zeros(1, 3);
            end

        case 'cauchy' % cauchy only
            peak_function = @cauchy; % Identify peaks as cauchy
            guess_params = zeros(max_n_peaks, 3);
            flat_spec = flat_iter;
            for guess = 1:max_n_peaks
                
                max_ind = find(flat_iter == max(flat_iter));
                max_height = flat_iter(max_ind);

                if max_height <= peak_threshold * std(flat_iter)
                    break
                end

                guess_freq = freqs(max_ind);
                guess_height = max_height;

                if guess_height <= min_peak_height
                    break
                end

                half_height = 0.5 * max_height;
                le_ind = sum(flat_iter(1:max_ind) <= half_height);
                ri_ind = length(flat_iter) - sum(flat_iter(max_ind:end) <= half_height);
                short_side = min(abs([le_ind,ri_ind]-max_ind));

                % Estimate gamma from FWHM. Calculate FWHM, converting to Hz, get guess gamma from FWHM
                fwhm = short_side * 2 * (freqs(2)-freqs(1));
                guess_gamma = fwhm/2;

                % Check that guess gamma isn't outside preset limits; restrict if so.
                % Note: without this, curve_fitting fails if given guess > or < bounds.
                if guess_gamma < gauss_std_limits(1)
                    guess_gamma = gauss_std_limits(1);
                end
                if guess_gamma > gauss_std_limits(2)
                    guess_gamma = gauss_std_limits(2);
                end

                % Collect guess parameters.
                guess_params(guess,:) = [guess_freq(1), guess_height, guess_gamma];

                % Subtract best-guess cauchy.
                peak_cauchy = cauchy(freqs, guess_freq(1), guess_height, guess_gamma);
                flat_iter = flat_iter - peak_cauchy;
            end

            guess_params(guess_params(:,1) == 0,:) = [];
            guess_params = drop_peak_cf(guess_params, proxThresh, [min(freqs) max(freqs)]);
            guess_params = drop_peak_overlap(guess_params, proxThresh);

            % If there are peak guesses, fit the peaks, and sort results.
            if ~isempty(guess_params)
                model_params = fit_peak_guess(guess_params, freqs, flat_spec, 2, guess_weight, gauss_std_limits,hOT);
            else
                model_params = zeros(1, 3);
            end

    end

end

%% Drop peak CF
function guess = drop_peak_cf(guess, bw_std_edge, freq_range)

    %% Notes
    %       Check whether to drop peaks based on center's proximity to the edge of the spectrum.
    %
    %       Parameters
    %       ----------
    %       guess : mx3 array, where m = No. of peaks.
    %           Guess parameters for peak fits.
    %
    %       Returns
    %       -------
    %       guess : qx3 where q <= m No. of peaks.
    %           Guess parameters for peak fits.

    %%
    cf_params = guess(:,1)';
    bw_params = guess(:,3)' * bw_std_edge;

    %% Check if peaks within drop threshold from the edge of the frequency range.
    keep_peak = abs(cf_params-freq_range(1)) > bw_params & ...
        abs(cf_params-freq_range(2)) > bw_params;

    %% Drop peaks that fail the center frequency edge criterion
    guess = guess(keep_peak,:);
end

%% Drop peak overlap
function guess = drop_peak_overlap(guess, proxThresh)

    %% Notes
    %       Checks whether to drop gaussians based on amount of overlap.
    %
    %       Parameters
    %       ----------
    %       guess : mx3 array, where m = No. of peaks.
    %           Guess parameters for peak fits.
    %       proxThresh: double
    %           Proximity threshold (in st. dev. or gamma) between two peaks.
    %
    %       Returns
    %       -------
    %       guess : qx3 where q <= m No. of peaks.
    %           Guess parameters for peak fits.
    %
    %       Note
    %       -----
    %       For any gaussians with an overlap that crosses the threshold,
    %       the lowest height guess guassian is dropped.

    %% Sort the peak guesses, so can check overlap of adjacent peaks
    guess = sortrows(guess);

    %% Calculate standard deviation bounds for checking amount of overlap
    bounds = [guess(:,1) - guess(:,3) * proxThresh, ...
        guess(:,1), guess(:,1) + guess(:,3) * proxThresh];

    %% Loop through peak bounds, comparing current bound to that of next peak
    drop_inds =  [];
    for ind = 1:size(bounds,1)-1
        b_0 = bounds(ind,:);
        b_1 = bounds(ind + 1,:);

        % Check if bound of current peak extends into next peak
        if b_0(2) > b_1(1)
            % If so, get the index of the gaussian with the lowest height (to drop)
            drop_inds = [drop_inds (ind - 1 + find(guess(ind:ind+1,2) == ...
                min(guess(ind,2),guess(ind+1,2))))];
        end
    end

    %% Drop any peaks guesses that overlap too much, based on threshold.
    guess(drop_inds,:) = [];
end

%% Core Models
function ys = gaussian(freqs, mu, hgt, sigma)

    %% Notes
    %       Gaussian function to use for fitting.
    %
    %       Parameters
    %       ----------
    %       freqs : 1xn array
    %           Frequency vector to create gaussian fit for.
    %       mu, hgt, sigma : doubles
    %           Parameters that define gaussian function (centre frequency,
    %           height, and standard deviation).
    %
    %       Returns
    %       -------
    %       ys :    1xn array
    %       Output values for gaussian function.
    ys = hgt*exp(-(((freqs-mu)./sigma).^2) /2);
end

function ys = cauchy(freqs, ctr, hgt, gam)

    %% Notes
    %       Cauchy function to use for fitting.
    %
    %       Parameters
    %       ----------
    %       freqs : 1xn array
    %           Frequency vector to create cauchy fit for.
    %       ctr, hgt, gam : doubles
    %           Parameters that define cauchy function (centre frequency,
    %           height, and "standard deviation" [gamma]).
    %
    %       Returns
    %       -------
    %       ys :    1xn array
    %       Output values for cauchy function.

    %%
    ys = hgt./(1+((freqs-ctr)/gam).^2);
end

% Exponential function for fitting
function ys = expo_function(freqs,params)

    %% Notes
    %       Exponential function to use for fitting 1/f, with a 'knee' (maximum at low frequencies).
    %
    %       Parameters
    %       ----------
    %       freqs : 1xn array
    %           Input x-axis values.
    %       params : 1x3 array (offset, knee, exp)
    %           Parameters (offset, knee, exp) that define Lorentzian function:
    %           y = 10^offset * (1/(knee + x^exp))
    %
    %       Returns
    %       -------
    %       ys :    1xn array
    %           Output values for exponential function.

    %%
    ys = params(1) - log10(abs(params(2)) +freqs.^params(3));
end

function ys = expo_nk_function(freqs, params)

    %% Notes
    %       Exponential function to use for fitting 1/f, without a 'knee'.
    %
    %       Parameters
    %       ----------
    %       freqs : 1xn array
    %           Input x-axis values.
    %       params : 1x2 array (offset, exp)
    %           Parameters (offset, exp) that define Lorentzian function:
    %           y = 10^offset * (1/(x^exp))
    %
    %       Returns
    %       -------
    %       ys :    1xn array
    %           Output values for exponential (no-knee) function.

    %%
    ys = params(1) - log10(freqs.^params(2));
end

%% Fit peak guess
function peak_params = fit_peak_guess(guess, freqs, flat_spec, peak_type, guess_weight, std_limits, hOT)

    %% Notes
        %     Fits a group of peak guesses with a fit function.
    %
    %     Parameters
    %     ----------
    %       guess : mx3 array, where m = No. of peaks.
    %           Guess parameters for peak fits.
    %       freqs : 1xn array
    %           Frequency values for the power spectrum, in linear scale.
    %       flat_iter : 1xn array
    %           Flattened (aperiodic removed) power spectrum.
    %       peakType : {'gaussian', 'cauchy', 'best'}
    %           Which types of peaks are being fitted.
    %       guess_weight : 'none', 'weak', 'strong'
    %           Parameter to weigh initial estimates during optimization.
    %       std_limits: 1x2 array
    %           Minimum and maximum standard deviations for distribution.
    %       hOT : 0 or 1
    %           Defines whether to use constrained optimization, fmincon, or
    %           basic simplex, fminsearch.
    %
    %       Returns
    %       -------
    %       peak_params : mx3, where m =  No. of peaks.
    %           Peak parameters post-optimization.
    
    %%
    if hOT % Use OptimToolbox for fmincon
        options = optimset('Display', 'off', 'TolX', 1e-3, 'TolFun', 1e-5, ...
            'MaxFunEvals', 3000, 'MaxIter', 3000); % Tuned options
        lb = [max([ones(size(guess,1),1).*freqs(1) guess(:,1)-guess(:,3)*2],[],2),zeros(size(guess(:,2))),ones(size(guess(:,3)))*std_limits(1)];
        ub = [min([ones(size(guess,1),1).*freqs(end) guess(:,1)+guess(:,3)*2],[],2),inf(size(guess(:,2))),ones(size(guess(:,3)))*std_limits(2)];
        peak_params = fmincon(@error_model_constr,guess,[],[],[],[], ...
            lb,ub,[],options,freqs,flat_spec, peak_type);
    else % Use basic simplex approach, fminsearch, with guess_weight
        options = optimset('Display', 'off', 'TolX', 1e-4, 'TolFun', 1e-5, ...
            'MaxFunEvals', 5000, 'MaxIter', 5000);
        peak_params = fminsearch(@error_model,...
            guess, options, freqs, flat_spec, peak_type, guess, guess_weight);
    end
end

%% Remove outliers
function SPRiNT = remove_outliers(SPRiNT,peak_function,opt)

    %% Notes
    % Helper function to remove outlier peaks according to user-defined specifications
    % Author: Luc Wilson

    %%
    timeRange = opt.maxtime.*opt.WinLength.*(1-opt.WinOverlap./100);
    nC = length(SPRiNT.channel);

    for c = 1:nC
        ts = [SPRiNT.channel(c).data.time];
        remove = 1;
        while any(remove)
            remove = zeros(length([SPRiNT.channel(c).peaks]),1);
            for p = 1:length([SPRiNT.channel(c).peaks])
                if sum((abs([SPRiNT.channel(c).peaks.time] - SPRiNT.channel(c).peaks(p).time) <= timeRange) &...
                        (abs([SPRiNT.channel(c).peaks.center_frequency] - SPRiNT.channel(c).peaks(p).center_frequency) <= opt.maxfreq)) < opt.minnear+1 % includes current peak
                    remove(p) = 1;
                end
            end
            SPRiNT.channel(c).peaks(logical(remove)) = [];
        end

        for t = 1:length(ts)

            if SPRiNT.channel(c).data(t).peak_params(1) == 0
                continue % never any peaks to begin with
            end
            p = [SPRiNT.channel(c).peaks.time] == ts(t);
            if sum(p) == size(SPRiNT.channel(c).data(t).peak_params,1)
                continue % number of peaks has not changed
            end
            peak_fit = zeros(size(SPRiNT.freqs));
            if any(p)
                SPRiNT.channel(c).data(t).peak_params = [[SPRiNT.channel(c).peaks(p).center_frequency]' [SPRiNT.channel(c).peaks(p).amplitude]' [SPRiNT.channel(c).peaks(p).st_dev]'];
                peak_pars = SPRiNT.channel(c).data(t).peak_params;
                for peak = 1:size(peak_pars,1)
                    peak_fit = peak_fit + peak_function(SPRiNT.freqs,peak_pars(peak,1),...
                        peak_pars(peak,2),peak_pars(peak,3));
                end
                ap_spec = log10(SPRiNT.channel(c).data(t).power_spectrum) - peak_fit;
                ap_pars = simple_ap_fit(SPRiNT.freqs, ap_spec, opt.aperiodic_mode, SPRiNT.channel(c).data(t).aperiodic_params(end));
                ap_fit = gen_aperiodic(SPRiNT.freqs, ap_pars, opt.aperiodic_mode);
                MSE = sum((ap_spec - ap_fit).^2)/length(SPRiNT.freqs);
                rsq_tmp = corrcoef(ap_spec+peak_fit,ap_fit+peak_fit).^2;
                % Return FOOOF results
                ap_pars(2) = abs(ap_pars(2));
                SPRiNT.channel(c).data(t).ap_fit = 10.^(ap_fit);
                SPRiNT.channel(c).data(t).fooofed_spectrum = 10.^(ap_fit+peak_fit);
                SPRiNT.channel(c).data(t).peak_fit = 10.^(peak_fit);
                SPRiNT.aperiodic_models(c,t,:) = SPRiNT.channel(c).data(t).ap_fit;
                SPRiNT.SPRiNT_models(c,t,:) = SPRiNT.channel(c).data(t).fooofed_spectrum;
                SPRiNT.peak_models(c,t,:) = SPRiNT.channel(c).data(t).peak_fit;
                SPRiNT.channel(c).aperiodics(t).offset = ap_pars(1);
                if length(ap_pars)>2 % Legacy FOOOF alters order of parameters
                    SPRiNT.channel(c).aperiodics(t).exponent = ap_pars(3);
                    SPRiNT.channel(c).aperiodics(t).knee_frequency = ap_pars(2);
                else
                    SPRiNT.channel(c).aperiodics(t).exponent = ap_pars(2);
                end
                SPRiNT.channel(c).stats(t).MSE = MSE;
                SPRiNT.channel(c).stats(t).r_squared = rsq_tmp(2);
                SPRiNT.channel(c).stats(t).frequency_wise_error = abs(ap_spec-ap_fit);

            else
                SPRiNT.channel(c).data(t).peak_params = [0 0 0];
                ap_spec = log10(SPRiNT.channel(c).data(t).power_spectrum) - peak_fit;
                ap_pars = simple_ap_fit(SPRiNT.freqs, ap_spec, opt.aperiodic_mode, SPRiNT.channel(c).data(t).aperiodic_params(end));
                ap_fit = gen_aperiodic(SPRiNT.freqs, ap_pars, opt.aperiodic_mode);
                MSE = sum((ap_spec - ap_fit).^2)/length(SPRiNT.freqs);
                rsq_tmp = corrcoef(ap_spec+peak_fit,ap_fit+peak_fit).^2;
                % Return FOOOF results
                ap_pars(2) = abs(ap_pars(2));
                SPRiNT.channel(c).data(t).ap_fit = 10.^(ap_fit);
                SPRiNT.channel(c).data(t).fooofed_spectrum = 10.^(ap_fit+peak_fit);
                SPRiNT.channel(c).data(t).peak_fit = 10.^(peak_fit);
                SPRiNT.aperiodic_models(c,t,:) = SPRiNT.channel(c).data(t).ap_fit;
                SPRiNT.SPRiNT_models(c,t,:) = SPRiNT.channel(c).data(t).fooofed_spectrum;
                SPRiNT.peak_models(c,t,:) = SPRiNT.channel(c).data(t).peak_fit;
                SPRiNT.channel(c).data(t).error = MSE;
                SPRiNT.channel(c).data(t).r_squared = rsq_tmp(2);
                SPRiNT.channel(c).aperiodics(t).offset = ap_pars(1);
                if length(ap_pars)>2 % Legacy FOOOF alters order of parameters
                    SPRiNT.channel(c).aperiodics(t).exponent = ap_pars(3);
                    SPRiNT.channel(c).aperiodics(t).knee_frequency = ap_pars(2);
                else
                    SPRiNT.channel(c).aperiodics(t).exponent = ap_pars(2);
                end
                SPRiNT.channel(c).stats(t).MSE = MSE;
                SPRiNT.channel(c).stats(t).r_squared = rsq_tmp(2);
                SPRiNT.channel(c).stats(t).frequency_wise_error = abs(ap_spec-ap_fit);
            end
        end
    end
end

%% Cluster peaks across time
function oS = cluster_peaks_dynamic2(oS)

    %% Notes
    % Second Helper function to cluster peaks across time
    % Author: Luc Wilson (2022)

    %%
    pthr = oS.options.proximity_threshold;
    for chan = 1:length(oS.channel)
        clustLead = [];
        nCl = 0;
        oS.channel(chan).clustered_peaks = struct();
        times = unique([oS.channel(chan).peaks.time]);
        all_peaks = oS.channel(chan).peaks;
        for time = 1:length(times)
            time_peaks = all_peaks([all_peaks.time] == times(time));
            % Initialize first clusters
            if time == 1
                nCl = length(time_peaks);
                for Cl = 1:nCl
                    oS.channel(chan).clustered_peaks(Cl).cluster = Cl;
                    oS.channel(chan).clustered_peaks(Cl).peaks(Cl) = time_peaks(Cl);
                    clustLead(Cl,1) = time_peaks(Cl).time;
                    clustLead(Cl,2) = time_peaks(Cl).center_frequency;
                    clustLead(Cl,3) = time_peaks(Cl).amplitude;
                    clustLead(Cl,4) = time_peaks(Cl).st_dev;
                    clustLead(Cl,5) = Cl;
                end
                continue
            end

            % Cluster "drafting stage": find points that make good matches to each cluster.
            for Cl = 1:nCl
                match = abs([time_peaks.center_frequency]-clustLead(Cl,2))./clustLead(Cl,4) < pthr;
                idx_tmp = find(match);
                if any(match)
                    % Add the best candidate peak
                    % Note: Auto-adds only peaks, but adds best candidate
                    % for multiple options
                    [tmp,idx] = min(([time_peaks(match).center_frequency] - clustLead(Cl,2)).^2 +...
                        ([time_peaks(match).amplitude] - clustLead(Cl,3)).^2 +...
                        ([time_peaks(match).st_dev] - clustLead(Cl,4)).^2);
                    oS.channel(chan).clustered_peaks(clustLead(Cl,5)).peaks(length(oS.channel(chan).clustered_peaks(clustLead(Cl,5)).peaks)+1) = time_peaks(idx_tmp(idx));
                    clustLead(Cl,1) = time_peaks(idx_tmp(idx)).time;
                    clustLead(Cl,2) = time_peaks(idx_tmp(idx)).center_frequency;
                    clustLead(Cl,3) = time_peaks(idx_tmp(idx)).amplitude;
                    clustLead(Cl,4) = time_peaks(idx_tmp(idx)).st_dev;
                    % Don't forget to remove the candidate from the pool
                    time_peaks(idx_tmp(idx)) = [];
                end
            end
            % Remaining peaks get sorted into their own clusters
            if ~isempty(time_peaks)
                for peak = 1:length(time_peaks)
                    nCl = nCl + 1;
                    Cl = nCl;
                    clustLead(Cl,1) = time_peaks(peak).time;
                    clustLead(Cl,2) = time_peaks(peak).center_frequency;
                    clustLead(Cl,3) = time_peaks(peak).amplitude;
                    clustLead(Cl,4) = time_peaks(peak).st_dev;
                    clustLead(Cl,5) = Cl;
                    oS.channel(chan).clustered_peaks(Cl).cluster = Cl;
                    oS.channel(chan).clustered_peaks(Cl).peaks(length(oS.channel(chan).clustered_peaks(clustLead(Cl,5)).peaks)+1) = time_peaks(peak);
                end
            end
            % Sort clusters based on most recent
            clustLead = sortrows(clustLead,1,'descend');
        end
    end
end

%% Error functions
function err = error_expo_function(params,xs,ys)
ym = expo_function(xs,params);
err = sum((ys - ym).^2);
end

function err = error_model_constr(params, xVals, yVals, peak_type)
    fitted_vals = 0;
    for set = 1:size(params,1)
        switch (peak_type)
            case 1 % Gaussian
                fitted_vals = fitted_vals + gaussian(xVals, params(set,1), params(set,2), params(set,3));
            case 2 % Cauchy
            fitted_vals = fitted_vals + cauchy(xVals, params(set,1), params(set,2), params(set,3));
        end
    end
    err = sum((yVals - fitted_vals).^2);
end