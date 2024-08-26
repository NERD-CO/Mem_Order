classdef MEM_ORDER < handle

    %% Properties
    properties(SetAccess = private)
        file
        subID
        info
        task
        LFP_data
        LFP_timestamps
        LFP_converted_timestamps
        LFP_idx
        epoch_idx
        fs
        freqs
        chanID
        chanIdx
        chanHemi
        chanSname
        brTABLE
        RegionLFP
        eventTimes
        eventIDs
        boundary
        boundaryAnalyzed
        boundaryIdx
        keepIdx
        wireID
        referencedData
        referencedChName
        crossStartIdx
        crossEndIdx
        clipOnIdx
        clipOffIdx
        questionIdx
        responseIdx
        fixation
        presentation
        presentation_corrected
        presentationPower
        question
        question_corrected
        questionPower
        ranksumprob
        fcMethod
        presentation_fc
        question_fc
        ch_names_corrected
        fcprob
    end

    methods

        %% Constructor
        function obj = MEM_ORDER(fname,task)

            %% Add path to NWB files
            addpath('C:\Users\Kevin_Tyner\Documents\MATLAB\matnwb-2.6.0.2');

            %% Make sure the file exists
            if isfile(fname)

                %% Add file name and task
                obj.file = fname;
                obj.task = task;

                %% Get subject ID
                obj.subID = extractBefore(extractAfter(fname,'DataFolder\'),'\NWB');

                %% Print
                fprintf('Behavioral data loaded for sub %s..\n',obj.subID)

                %% Get event data
                folder = extractBefore(fname,'NWBProcessing');
                task_info = append(folder,'Behavioral_Data\Raw\');
                files = dir(task_info);
                files = {files(~[files.isdir]).name}';
                I = contains(files,'.mat') & contains(files,task);
                file = string(append(task_info,files(I)));
                if isfile(file)
                    obj.info = load(file);
                end

                %% Load the data
                data = nwbRead(fname);

                %% Load macrowire timestamps
                timestamps = data.processing.get('ecephys').nwbdatainterface.get...
                    ('LFP').electricalseries.get('MacroWireSeries').timestamps.load;

                %% Downsample timestamps
                obj.LFP_timestamps = downsample(timestamps,8);

                %% Convert timestamps
                obj.LFP_converted_timestamps = (obj.LFP_timestamps - obj.LFP_timestamps(1,1))/(1e6); % convert from microsec to sec

                %% Get voltage data for all macrowires and their time stamps
                temp_LFP_data = double(data.processing.get('ecephys').nwbdatainterface.get...
                    ('LFP').electricalseries.get('MacroWireSeries').data.load); % do the conversion before adding to object

                %% Get the sampling frequency
                LFP_sessionInfo = data.processing.get('ecephys').nwbdatainterface.get...
                    ('LFP').electricalseries.get('MacroWireSeries');
                obj.fs = str2double(cell2mat(extractBetween(LFP_sessionInfo.description,'= ',':')));

                %% Convert LFP data and add to object
                obj.LFP_data = temp_LFP_data .* LFP_sessionInfo.data_conversion;

                %% Channel IDs
                chanLabels = cellstr(data.general_extracellular_ephys_electrodes.vectordata.get('label').data.load()); %use MA only!
                MAchan = find(contains(chanLabels,'MA_'));
                chanID = cellstr(data.general_extracellular_ephys_electrodes.vectordata.get('location').data.load());
                hemisphere = cellstr(data.general_extracellular_ephys_electrodes.vectordata.get('hemisph').data.load());
                shortBnames = cellstr(data.general_extracellular_ephys_electrodes.vectordata.get('shortBAn').data.load());
                wireID = data.general_extracellular_ephys_electrodes.vectordata.get('channID').data.load();

                obj.chanID = chanID(MAchan);
                obj.chanHemi = hemisphere(MAchan);
                obj.chanSname = shortBnames(MAchan);
                obj.wireID = wireID(MAchan);

                %% Find unique brain regions and turn into a table
                tmpBRegUni = unique(obj.chanSname);
                hemiTemp = cell(length(tmpBRegUni),1);
                longName = cell(length(tmpBRegUni),1);
                for ui = 1:length(tmpBRegUni)
                    tmpIND = find(matches(obj.chanSname,tmpBRegUni{ui}),1,'first');
                    hemiTemp{ui} = obj.chanHemi{tmpIND};
                    longName{ui} = obj.chanID{tmpIND};
                end
                obj.brTABLE = table(tmpBRegUni,hemiTemp,longName,'VariableNames',{'SEEGele',...
                    'Hemisphere','LongBRname'});

                %% Assign LFP to each brain region
                obj.RegionLFP = cell(length(obj.brTABLE.SEEGele),1);
                for k = 1:height(obj.brTABLE)
                    reg = strcmp(obj.brTABLE.SEEGele{k},obj.chanSname);
                    obj.RegionLFP{k} = obj.LFP_data(reg,:);
                end

                %% Extract event key
                obj.eventTimes = data.acquisition.get('events').timestamps.load();
                temp_eventIDs = cellstr(data.acquisition.get('events').data.load());

                %% Convert eventIDs from hexadecimal
                I = contains(temp_eventIDs,'TTL');
                obj.eventTimes = obj.eventTimes(I);
                temp_eventIDs2 = temp_eventIDs(I);
                TTL_task = extractBetween(temp_eventIDs2,'(',')');
                obj.eventIDs = cellfun(@(x) hex2dec(x),TTL_task,'UniformOutput',true);

            end
        end

        %% Methods
        
            %% Grab channels of interest
            function selectChannels(obj,channels,hemi)

                %% Notes
                % This function selects the data that corresponds to
                % channels of interest for the analysis.  At the cuttent
                % iteration, this function only handels one brain region at
                % a time.  Inputs for channels are strings such as
                % "anterior hippocampus" or "amygdala", and hemi are "L" or "R".

                %% Find channels of interest
                I = strcmp(obj.chanID,channels) & strcmp(obj.chanHemi,hemi);

                %% Remove channels
                obj.LFP_data = obj.LFP_data(I,:);
                obj.chanID = obj.chanID(I,:);
                obj.chanHemi = obj.chanHemi(I,:);
                obj.chanSname = obj.chanSname(I,:);
                obj.wireID = obj.wireID(I,:);

                %% Modify channel names
                for z = 1:length(obj.chanID)
                    obj.chanIdx{z,1} = strcat(obj.chanSname{z,1},'_',num2str(z));
                end
            end

            %% Reject bad channels
            function identifyBads(obj)

                %% Notes
                % This function checks if channels are by by determining if
                % 20% or more of the data in the channel is greater than or
                % less than the mean of each channel +/- 2 standard
                % deviations.

                %% Loop through channels to identify bads
                for ii = 1:length(obj.chanIdx)

                    %% Calculate mean and standard deviation
                    ch_mean = mean(obj.LFP_data(ii,:));
                    ch_std = std(obj.LFP_data(ii,:));

                    %% Set thresholds
                    upper_th = ch_mean + 2*ch_std;
                    lower_th = ch_mean - 2*ch_std;

                    %% Determine if channel is by by comparing to thresholds
                    if sum(obj.LFP_data(ii,:) > upper_th | obj.LFP_data(ii,:) < lower_th) > 0.2*length(obj.LFP_data(ii,:))
                        obj.keepIdx(ii,1) = 0;
                    else
                        obj.keepIdx(ii,1) = 1;
                    end
                end

                %% Convert to logical
                obj.keepIdx = logical(obj.keepIdx);

                %% Replace bad channels with NaN
                obj.LFP_data(~obj.keepIdx,:) = NaN;

            end

            %% Perform bipolar referencing
            function bipolarMontage(obj)

                %% Notes
                % This function performs bipolar montage referencing on the
                % channels that made it through data selection.  In the
                % event a channel id dropped and the remaining channels are
                % non-consecutive (e.g. A_1, A_2, A_4), the method will
                % reference to the next available channel.

                %% Loop through and grab the data associated with good channels
                ii = 1;
                temp = cell(length(obj.keepIdx),1);
                ch_name = cell(length(obj.keepIdx),1);
                while ii <= length(obj.keepIdx)
                    if obj.keepIdx(ii,1) == 1
                        temp{ii,1} = obj.LFP_data(ii,:);
                        ch_name{ii,1} = obj.chanIdx{ii,1};
                        ii = ii + 1;
                    else
                        next_idx = find(obj.keepIdx(ii+1,1) == 1,1,"first");
                        if ~isempty(next_idx)
                            ii = ii + next_idx;
                            temp{ii,1} = obj.LFP_data(ii,:);
                            ch_name{ii,1} = obj.chanIdx{ii,1};
                        else
                            break
                        end
                    end
                end

                %% Find empty cells and fix wire IDs
                I = cellfun(@isempty,temp);
                obj.wireID = obj.wireID(~I,:);

                %% Remove empty cells
                temp = temp(~cellfun('isempty',temp));
                ch_name = ch_name(~cellfun('isempty',ch_name));

                %% Convert
                temp = cell2mat(temp);

                %% Loop to perform bipolar referencing
                for ii = 1:length(ch_name) - 1
                    obj.referencedData(ii,:) = temp(ii,:) - temp(ii+1,:);
                    obj.referencedChName{ii,1} = strcat(ch_name{ii,1},'-',ch_name{ii+1,1});
                end

            end

            %% Identify TTLs
            function identifyTTLs(obj)

                %% Notes

                %% Remove bad TTLs
                I = obj.eventIDs == 11 | obj.eventIDs == 1 | obj.eventIDs == 2 | obj.eventIDs == 3;
                obj.eventIDs = obj.eventIDs(I);
                obj.eventTimes = obj.eventTimes(I);

                %% Find LFP idx
                for ii = 1:length(obj.eventIDs)
                    [~,obj.LFP_idx(ii,1)] = min(abs(obj.eventTimes(ii,1) - obj.LFP_timestamps));
                end

                %% Preallocate
                obj.epoch_idx = zeros(length(obj.eventIDs),1);
                start = 1;

                %% Assign to an epoch
                for jj = 1:length(obj.eventIDs)
                    if (obj.eventIDs(ii,1) == 11) && (obj.eventIDs(ii+1,1) == 1) && (obj.eventIDs(ii+2,1) == 2) && ...
                            (obj.eventIDs(ii+3,1) == 3)
                        obj.epoch_idx(jj,1) = start;
                        start = start + 1;
                    end
                end

                %% Get the field name
                fieldNames = fieldnames(obj.info);
                for z = 1:length(fieldNames)
                    if contains(fieldNames(z,1),'respMat')
                        matchedField = fieldNames(z,1);
                    end
                end

                %% Calculate time differentials
                crossDiff = zeros(length(obj.info.(matchedField{1,1})),1);
                quesDiff = zeros(length(obj.info.(matchedField{1,1})),1);
                for ii = 1:length(obj.info.(matchedField{1,1}))
                    crossDiff(ii,1) = obj.info.(matchedField{1,1})(ii).CrossEnd - obj.info.(matchedField{1,1})(ii).CrossStart;
                    quesDiff(ii,1) = obj.info.(matchedField{1,1})(ii).respTime - obj.info.(matchedField{1,1})(ii).QuesStart;
                end

                %% Convert from time to samples
                crossDiff = floor(crossDiff .* obj.fs);
                quesDiff = floor(quesDiff .* obj.fs);

                %% Move indices to the object
                obj.crossStartIdx = obj.LFP_idx(obj.eventIDs == 11);
                obj.crossEndIdx = obj.crossStartIdx + crossDiff;
                obj.clipOnIdx = obj.LFP_idx(obj.eventIDs == 1);
                obj.clipOffIdx = obj.LFP_idx(obj.eventIDs == 2);
                obj.questionIdx = obj.LFP_idx(obj.eventIDs == 3);
                obj.responseIdx = obj.questionIdx + quesDiff;

            end

            %% Identify task boundaries
            function identifyBoundaries(obj)

                %% Notes
                % This function will determine the boundary type for each
                % trial in the task.

                %% Get the field name
                fieldNames = fieldnames(obj.info);
                for z = 1:length(fieldNames)
                    if contains(fieldNames(z,1),'respMat')
                        matchedField = fieldNames(z,1);
                    end
                end

                %% Find appropriate field name for clip/frame
                x = fieldnames(obj.info.(matchedField{1,1}));
                I = strcmp(x,'FrameName') | strcmp(x,'ClipName');

                %% Loop and identify each boundary
                obj.boundary = cell(length(obj.info.(matchedField{1,1})),1);
                for ii = 1:length(obj.info.(matchedField{1,1}))
                    obj.boundary{ii,1} = extractBefore(obj.info.(matchedField{1,1})(ii).(x{I}),'_');
                end
            end

            %% Grab the data
            function allocateData(obj)

                %% Notes
                % This function separates the data into the three phases;
                % fixation, presentation, and question, and adds that data
                % to the object.

                %% Create anonymous function handle
                extractSegments = @(s,e) obj.referencedData(:,s:e-1);

                %% Grab the data
                obj.fixation = arrayfun(@(i) extractSegments(obj.crossStartIdx(i),obj.crossEndIdx(i)), 1:length(obj.boundary), 'UniformOutput', false )';
                obj.presentation = arrayfun(@(i) extractSegments(obj.clipOnIdx(i),obj.clipOffIdx(i)), 1:length(obj.boundary), 'UniformOutput', false )';
                obj.question = arrayfun(@(i) extractSegments(obj.questionIdx(i),obj.responseIdx(i)), 1:length(obj.boundary), 'UniformOutput', false )';
                
            end

            %% Baseline subtraction
            function baselineCorrection(obj)

                %% Notes
                % This function calculates the mean amplitude for each
                % channel during the fixation phase and subtracts it from
                % each value for every channel in the presentation and
                % question phase.

                %% Preallocate
                obj.presentation_corrected = cell(length(obj.presentation),1);
                obj.question_corrected = cell(length(obj.question),1);

                %% Loop for baseline correction
                for ii = 1:length(obj.fixation)

                    %% Calculate channel means
                    ch_means = mean(obj.fixation{ii,1},2);

                    %% Subtract the mean from presentation data
                    if ~isempty(obj.presentation{ii,1})
                        obj.presentation_corrected{ii,1} = obj.presentation{ii,1} - ch_means;
                    end

                    %% Subtract the mean from question data
                    if ~isempty(obj.question{ii,1})
                        obj.question_corrected{ii,1} = obj.question{ii,1} - ch_means;
                    end
                end
            end

            %% Plot regions
            function plotPhases(obj)

                %% Notes
                % This function plots the baseline corrected LFP and shades
                % the regions that correspond to the fixation, presentation
                % and question phases.

                %% Plot
                figure;
                hold on
                for z = 1:size(obj.referencedData,1)
                    plot(obj.referencedData(z,:))
                end

                %% Bands for cross presentation
                bands = [obj.crossStartIdx,obj.crossEndIdx];
                xp = [bands fliplr(bands)];
                yp = ([[1;1]*min(ylim); [1;1]*max(ylim)]*ones(1,size(bands,1))).';
                for k = 1:size(bands,1)
                    patch(xp(k,:),yp(k,:),[1 0 0],'FaceAlpha',0.3,'EdgeColor','none')
                end

                %% Bands for clip onset
                bands = [obj.clipOnIdx,obj.clipOffIdx];
                xp = [bands fliplr(bands)];
                for k = 1:size(bands,1)
                    patch(xp(k,:),yp(k,:),[0 1 0],'FaceAlpha',0.3,'EdgeColor','none')
                end

                %% Bands for question onset
                bands = [obj.questionIdx,obj.responseIdx];
                xp = [bands fliplr(bands)];
                for k = 1:size(bands,1)
                    patch(xp(k,:),yp(k,:),[0 0 1],'FaceAlpha',0.3,'EdgeColor','none')
                end

            end

            %% Compute power
            function computePower(obj,freq)

                %% Notes
                % This function computes the power on the presentation and
                % question phases of each task and ignores trials with
                % insufficient data.

                %% Add freq to object
                obj.freqs = freq;

                %% Loop through presentation
                for ii = 1:length(obj.presentation_corrected)

                    %% Check if there is data present
                    if ~isempty(obj.presentation_corrected{ii,1})

                        %% Grab the data
                        temp = obj.presentation_corrected{ii,1};

                        %% Loop through channels and compute CWT
                        for jj = 1:size(temp,1)

                            %% Check for NaNs
                            if ~any(isnan(temp(jj,1)))

                                %% Compute power
                                [wt,f] = cwt(temp(jj,:),obj.fs);

                                %% Identify power in range of interest
                                I = f >= freq(1,1) & f <= freq(2,1);
                                obj.presentationPower{ii,1}(jj,:) = mean(abs(wt(I,:)),1);
                            else
                                obj.presentationPower{ii,1}(jj,1:length(temp)) = NaN;
                            end
                        end
                    end
                end

                %% Loop through question
                for ii = 1:length(obj.question_corrected)

                    %% Check if there is data present
                    if ~isempty(obj.question_corrected{ii,1})

                        %% Grab the data
                        temp = obj.question_corrected{ii,1};

                        %% Loop through channels and compute CWT
                        for jj = 1:size(temp,1)

                            %% Check for NaNs
                            if ~any(isnan(temp(jj,1))) && size(temp,2) > obj.fs/5 % more than 200 ms of data

                                %% Compute power
                                [wt,f] = cwt(temp(jj,:),obj.fs);

                                %% Identify power in range of interest
                                I = f >= freq(1,1) & f <= freq(2,1);
                                if any(I)
                                    obj.questionPower{ii,1}(jj,:) = mean(abs(wt(I,:)),1);
                                else
                                    obj.questionPower{ii,1} = [];
                                end
                            else
                                %obj.questionPower{ii,1}(jj,1:size(temp,2)) = NaN;
                                obj.questionPower{ii,1} = [];
                            end
                        end
                    end
                end
            end

            %% Statistical comparison
            function computeStats(obj,boundary)

                %% Notes
                % Use first 200 ms of question phase to examine changes in power between
                % presentation and question.  The boundary argument uses
                % HB, SB, or NB and grabs only the trials that correspond
                % to the provided boundary.

                %% Find number of samples for 200 ms
                nSamp = obj.fs/5; % samples per 1000 ms - divide by 5 for samples per 200 ms

                %% Add boundary to object
                obj.boundaryAnalyzed = boundary;

                %% Find only the requested boundary
                I = zeros(length(obj.boundary),1);
                for ii = 1:length(obj.boundary)
                    x = strcmp(obj.boundary{ii,1},boundary);
                    if any(x)
                        I(ii) = 1;
                    else
                        I(ii) = 0;
                    end
                end
                I = logical(I);
                obj.presentationPower = obj.presentationPower(I);
                obj.questionPower = obj.questionPower(I);
                obj.boundaryIdx = I;

                %% Average power over the desired samples for presentation task
                avg_presentationPower = zeros(size(obj.referencedData,1),length(obj.presentationPower));
                for ii = 1:length(obj.presentationPower)
                    avg_presentationPower(:,ii) = mean(obj.presentationPower{ii,1},2);
                end

                %% Average power over the desired sample for the question task
                start = 1;
                avg_questionPower = zeros(size(obj.referencedData,1),length(obj.questionPower));
                for ii = 1:length(obj.questionPower)
                    if ~isempty(obj.questionPower{ii,1})
                        avg_questionPower(:,start) = mean(obj.questionPower{ii,1}(:,1:nSamp),2);
                        start = start + 1;
                    end
                end

                %% Do the power comparison
                for jj = 1:size(avg_presentationPower,1)
                    obj.ranksumprob(jj,1) = ranksum(avg_presentationPower(jj,:),avg_questionPower(jj,:));
                end

            end

            %% Functional Connectivity
            function functionalConnectivity(obj,method)

                %% Notes
                % This function calculates the function connectivity
                % between channels in the corrected presentation and
                % question data.  Acceptable methods are "granger" for
                % granger causality and "coherence" for magnitude-squared
                % coherence.

                %% Add to object
                obj.fcMethod = method;

                %% Loop to compute presentation functional connectivity
                for ii = 1:length(obj.presentation_corrected)
                    if ~isempty(obj.presentation_corrected{ii,1})

                        %% Get the data
                        temp = obj.presentation_corrected{ii,1};

                        %% Preallocate
                        obj.presentation_fc{ii,1} = zeros(size(temp,1),size(temp,1));

                        %% Check FC method
                        if strcmp(method,'granger')

                            %% Calculate number of lags needed
                            lags = ceil(1/min(obj.freqs) * obj.fs);

                            %% Exit if not enough data
                            if size(temp,2) < 4*lags
                                obj.presentation_fc{ii,1} = [];
                                continue
                            end

                            %% Loop over channels
                            for jj = 1:size(temp,1)
                                for kk = 1:size(temp,1)
                                    if jj ~= kk

                                        %% Grab channel data
                                        sig_reduced = temp(jj,:);
                                        sig_full = [temp(jj,:);temp(kk,:)];

                                        %% Calculate residuals
                                        [R_reduced,~] = AutoregressiveProcess(sig_reduced,lags);
                                        [R_full,~] = AutoregressiveProcess(sig_full,lags);

                                        %% Calculate variances
                                        var_reduced = var(R_reduced(1,:),0,2);
                                        var_full = var(R_full(1,:),0,2);

                                        %% Write FC value
                                        obj.presentation_fc{ii,1}(jj,kk) = log(var_reduced/var_full);

                                    end
                                end
                            end

                        elseif strcmp(method,'coherence')

                            %% Exit if not enough data
                            if size(temp,2) < 9
                                obj.presentation_fc{ii,1} = [];
                                continue
                            end

                            %% Compute mscohere
                            for jj = 1:size(temp,1)
                                for kk = 1:size(temp,1)
                                    if jj ~= kk

                                        %% Calculate magnitude squared coherence
                                        [x,f] = mscohere(temp(jj,:),temp(kk,:),[],[],[],obj.fs);

                                        %% Grab values in frequency range
                                        J = f >= min(obj.freqs) & f <= max(obj.freqs);

                                        %% Write the mean values
                                        obj.presentation_fc{ii,1}(jj,kk) = mean(x(J));
                                    elseif jj == kk
                                        obj.presentation_fc{ii,1}(jj,kk) = 0;
                                    end
                                end
                            end
                        end
                    end
                end
                
                %% Loop to compute question functional connectivity
                for ii = 1:length(obj.question_corrected)
                    if ~isempty(obj.question_corrected{ii,1})

                        %% Get the data
                        temp = obj.question_corrected{ii,1};

                        %% Preallocate
                        obj.question_fc{ii,1} = zeros(size(temp,1),size(temp,1));

                        %% Check FC method
                        if strcmp(method,'granger')

                            %% Calculate number of lags needed
                            lags = ceil(1/min(obj.freqs) * obj.fs);

                            %% Exit if not enough data
                            if size(temp,2) < 4*lags
                                obj.question_fc{ii,1} = [];
                                continue
                            end

                            %% Loop over channels
                            for jj = 1:size(temp,1)
                                for kk = 1:size(temp,1)
                                    if jj ~= kk

                                        %% Grab channel data
                                        sig_reduced = temp(jj,:);
                                        sig_full = [temp(jj,:);temp(kk,:)];

                                        %% Calculate residuals
                                        [R_reduced,~] = AutoregressiveProcess(sig_reduced,lags);
                                        [R_full,~] = AutoregressiveProcess(sig_full,lags);

                                        %% Calculate variances
                                        var_reduced = var(R_reduced(1,:),0,2);
                                        var_full = var(R_full(1,:),0,2);

                                        %% Write FC value
                                        obj.question_fc{ii,1}(jj,kk) = log(var_reduced/var_full);

                                    end
                                end
                            end
                        elseif strcmp(method,'coherence')

                            %% Exit if not enough data
                            if size(temp,2) < 9
                                obj.question_fc{ii,1} = [];
                                continue
                            end

                            %% Compute mscohere
                            for jj = 1:size(temp,1)
                                for kk = 1:size(temp,1)
                                    if jj ~= kk

                                        %% Calculate magnitude squared coherence
                                        [x,f] = mscohere(temp(jj,:),temp(kk,:),[],[],[],obj.fs);

                                        %% Grab values in frequency range
                                        J = f >= min(obj.freqs) & f <= max(obj.freqs);

                                        %% Write the mean values
                                        obj.question_fc{ii,1}(jj,kk) = mean(x(J));
                                    elseif jj == kk
                                        obj.question_fc{ii,1}(jj,kk) = 0;
                                    end
                                end
                            end
                        end
                    end
                end
            end

            %% Analyze functional connectivity
            function analyzeFC(obj)

                %% Notes
                % This function calculates the functional connectivity
                % during the presentation and question phase, and
                % identifies what channel pairs are significant after
                % Benjamini-Hochberg correction for multiple comparisons.
                % This function also produces figures for the mean
                % functional connectivity between the two phases.

                %% Grab FC matrices of interest
                obj.presentation_fc = obj.presentation_fc(obj.boundaryIdx);
                obj.question_fc = obj.question_fc(obj.boundaryIdx);

                %% Calculate the mean FC for each phase
                pres_fc = zeros(length(obj.referencedChName),length(obj.referencedChName));
                pres_count = 0;
                for ii = 1:length(obj.presentation_fc)
                    if ~isempty(obj.presentation_fc{ii,1})
                        pres_count = pres_count + 1;
                        pres_fc = pres_fc + obj.presentation_fc{ii,1};
                    else
                        continue
                    end
                end
                pres_fc = pres_fc./pres_count;

                ques_fc = zeros(length(obj.referencedChName),length(obj.referencedChName));
                ques_count = 0;
                for ii = 1:length(obj.question_fc)
                    if ~isempty(obj.question_fc{ii,1})
                        ques_count = ques_count + 1;
                        ques_fc = ques_fc + obj.question_fc{ii,1};
                    else
                        continue
                    end
                end
                ques_fc = ques_fc./ques_count;

                fc_presentation_max = max(max(pres_fc));
                fc_question_max = max(max(ques_fc));
                y = max(fc_question_max,fc_presentation_max);

                fc_presentation_min = min(min(pres_fc));
                fc_question_min = min(min(ques_fc));
                x = min(fc_presentation_min,fc_question_min);



                figure
                imagesc(pres_fc);
                clim manual;
                clim([x y]);
                cb = colorbar;
                cb.Label.String = 'Presentation FC';

                figure
                imagesc(ques_fc);
                clim manual;
                clim([x y])
                cb = colorbar;
                cb.Label.String = 'Question FC';

                %% Calculate probability values for FC matrices
                p = ones(length(obj.referencedChName)*length(obj.referencedChName),1);
                ch_names = cell(length(obj.referencedChName)*length(obj.referencedChName),1);
                ch_ref = cell(length(obj.referencedChName),length(obj.referencedChName));
                count = 1;
                for ii = 1:length(obj.referencedChName)
                    for jj = 1:length(obj.referencedChName)
                        pres_val = NaN(length(obj.question_fc),1);
                        ques_val = NaN(length(obj.question_fc),1);
                        ch_ref{ii,jj} = append(obj.referencedChName{ii,1},'/',obj.referencedChName{jj,1});
                        if ii ~= jj
                            for kk = 1:length(obj.question_fc)
                                if ~isempty(obj.question_fc{kk,1})
                                    pres_val(kk,1) = obj.presentation_fc{kk,1}(ii,jj);
                                    ques_val(kk,1) = obj.question_fc{kk,1}(ii,jj);
                                else
                                    continue
                                end
                            end
                            [~,p(count,1)] = ttest(pres_val,ques_val);
                            ch_names{count,1} = append(obj.referencedChName{ii,1},'/',obj.referencedChName{jj,1});
                            count = count + 1;
                        else
                            ch_names{count,1} = append(obj.referencedChName{ii,1},'/',obj.referencedChName{jj,1});
                            count = count + 1;
                            continue
                        end
                    end
                end

                %% Correct for False discovery rate with Benjamini-Hochberg
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
                obj.fcprob = val(corrected);
                obj.ch_names_corrected = new_names_corrected(idx);

            end

            %% Next function

    end
end