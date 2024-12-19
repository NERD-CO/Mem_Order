function [] = checkSACCADE_MO(gazeFILEn)

load(gazeFILEn, 'outInfo')

eyeIDs = {'leftEYE','rightEYE'};

for eeii = 1:length(eyeIDs) % Loop through var(i) fields

    % CLEAN up duplicates
    outInfo.(eyeIDs{eeii}).GAZEcl =...
        cleanUPdupSacs(outInfo.(eyeIDs{eeii}).GAZEcl);

    % CHECK SACCADES
    % 1. Sanity plot loop
    [newSACCADEcol] = plotSACCADEoverlay(outInfo.(eyeIDs{eeii}).GAZEcl,...
        outInfo.(eyeIDs{eeii}).gaze_Scr_time,...
        outInfo.(eyeIDs{eeii}).gaze_Scr);

    outInfo.(eyeIDs{eeii}).GAZEcl2 = transpose(newSACCADEcol);

    disp(['Eye ', num2str(eeii) ,' done'])
end % VARIANT field loop

save(gazeFILEn, 'outInfo')


end





% CLEAN UP DUPLICATE SACCADES

function [saccadeCLEAN] = cleanUPdupSacs(saccadeOLD)

saccadeCLEAN = saccadeOLD;
for trii = 1:height(saccadeOLD)
    tmpTRIAL = saccadeOLD{trii}.saccades;

    if ~istable(tmpTRIAL)
        continue
    end

    tabSTARTtms = tabulate(categorical(tmpTRIAL.starttime));

    % Find rows with greater than 1
    countROWS = [tabSTARTtms{:,2}];

    % First check for any duplicates
    if any(countROWS ~= 1)

        indexOfDups = countROWS ~= 1;
        numOfDups = sum(indexOfDups);

        for ndi = 1:numOfDups
            [tmpTRIAL] = getNEWtable(tmpTRIAL);
        end
        % UPDATE TRIAL TABLE
        saccadeCLEAN{trii}.saccades = tmpTRIAL;
    else
        continue
    end
end



end




function [newTable] = getNEWtable(oldTable)

tabSTARTtms = tabulate(categorical(oldTable.starttime));

% Find rows with greater than 1
countROWS = [tabSTARTtms{:,2}];

rowID = find(countROWS > 1,1,'first');
timeID = tabSTARTtms{rowID,1};

oldTableLOCS = find(ismember(oldTable.starttime,str2double(timeID)));

removeLOCS = oldTableLOCS(~ismember(oldTableLOCS,min(oldTableLOCS)));

newTable = oldTable;

newTable(removeLOCS,:) = [];


end


%%%% SACCADE CHECK

function [newSACCADEcolumn] = plotSACCADEoverlay(gazeCOLUMN,gazeTM,gazePT)

cmap2use = colormap("plasma");
close all

for i = 1:height(gazeCOLUMN)

    tmpTRIAL = gazeCOLUMN{i}.saccades;

    if ~istable(tmpTRIAL)
        newSACCADEcolumn{i}.saccades = nan;
        newSACCADEcolumn{i}.fixations = nan;
        continue
    else

        % Plot raw Scr points
        gazeRAW = gazePT{i};
        gazeTime = gazeTM{i};
        rawS = scatter(gazeRAW(:,1),gazeRAW(:,2),60,[0.5 0.5 0.5],'filled');
        rawS.MarkerEdgeColor = "none";
        rawS.MarkerFaceAlpha = 0.7;
        % Loop through saccade table

        params.velocityThreshold = 3; % multiplier relative to median velocity (adjust as needed)
        params.minDuration = 7; % Minimum duration in milliseconds (adjust as needed)
        params.minAmplitude = 0.5;  % Minimum saccade amplitude in pixels (adjust as needed)
        % params.lossWindow = 100; %ms
        % params.lossThreshold = 2; %variance of 2 pixels^2
        params.windowSize = 100;
        params.velocityThresholdMultiplier = 4;

        % sacPARMS - amplitude, duration, velocity, peakVelocity,
        % direction, starttime, endtime

        [saccades , rmsDATA] = detectSaccadesORG(gazeRAW, 1000, params);

        saccadeMetrics = computeSaccadeMetrics(gazeRAW, saccades, 1000);

        disp(num2str(i))

        if isempty(saccades)
            newSACCADEcolumn{i}.saccades = nan;
            newSACCADEcolumn{i}.fixations = nan;
            continue
        end


        starttime = gazeTime(saccades(:,1));
        endtime = gazeTime(saccades(:,2));

        newSACCADEtab = table(saccadeMetrics.amplitude, transpose(saccadeMetrics.duration),...
            saccadeMetrics.meanVelocity, saccadeMetrics.peakVelocity,...
            saccadeMetrics.direction, starttime , endtime, 'VariableNames',...
            {'amplitude','duration','velocity','peakVelocity','direction',...
            'starttime','endtime'});

        newSACCADEcolumn{i}.saccades = newSACCADEtab;
        newSACCADEcolumn{i}.fixations = nan;

        cmapSize = size(cmap2use, 1);
        indices = round(linspace(1, cmapSize, height(saccades))); % Generate evenly spaced indices
        spacedColors = cmap2use(indices, :);

        for ssi = 1:height(saccades)%1:height(tmpTRIAL)

            % if tmpTRIAL.duration(ssi) > 200
            %     continue
            % else

            hold on
            % stTIME = tmpTRIAL.starttime(ssi);
            % edTIME = tmpTRIAL.endtime(ssi);
            %
            % gazeTIME = gazeTM{i};
            % sacINDEX = [find(gazeTIME == stTIME),find(gazeTIME == edTIME)];
            % saccadePTS = gazeRAW(sacINDEX(1):sacINDEX(2),:);
            sacS = scatter(gazeRAW(saccades(ssi,1):saccades(ssi,2),1),gazeRAW(saccades(ssi,1):saccades(ssi,2),2),50,spacedColors(ssi,:));
            % sacS = scatter(saccadePTS(:,1),saccadePTS(:,2),50,spacedColors(ssi,:));
            % sacS.MarkerEdgeColor = "none";
            % sacS.MarkerFaceAlpha = 0.5;
            sacS.LineWidth = 1;

            % end
        end % SACCADE loop
    end % Make sure not NAN
    title(num2str(i))
    pause(0.01)
    close all
end % end of gaze column loop
end

% function [saccades, RMSdist] = detectSaccades_Dynamic(eyePositions, samplingRate, params)
% % DETECTSACCADES Detects saccades in eye position data using a dynamic threshold.
% %
% %   Inputs:
% %       eyePositions: 2D array (n x 2) of eye positions (x, y).
% %       samplingRate: Sampling rate of the data (Hz).
% %       params: Structure containing detection parameters:
% %           .velocityThresholdMultiplier: Velocity threshold multiplier relative to median velocity
% %           .minDuration: Minimum saccade duration (ms)
% %           .minAmplitude: Minimum saccade amplitude (pixels)
% %           .windowSize: Sliding window size (ms) for the median calculation.
% %
% %   Outputs:
% %       saccades: A matrix where each row is a saccade, with columns:
% %              start index, end index.
%
%
%     % Parameter handling (set defaults if no params passed)
%     if nargin < 3 || isempty(params)
%         params.velocityThresholdMultiplier = 4;
%         params.minDuration = 20; %ms
%         params.minAmplitude = 0.5; %pixels
%         params.windowSize = 100; % ms (window size)
%         warning('Parameter structure not found or empty; setting default saccade detection parameters.');
%     end
%
%
%     windowSizeSamples = round(params.windowSize / 1000 * samplingRate); %convert ms to samples
%
%     % 1. Calculate Velocities
%     velocities = calculateVelocity(eyePositions, samplingRate);
%     speeds = sqrt(sum(velocities.^2, 2)); % Magnitude of the velocities
%
%
%     % 2. Dynamic Thresholding (using a sliding window)
%     numPoints = length(speeds);
%     saccadeCandidates = zeros(size(speeds)); % start as zeros
%
%
%     for i = 1:numPoints
%         windowStart = max(1, i - floor(windowSizeSamples/2));
%         windowEnd = min(numPoints, i + floor(windowSizeSamples/2));
%         windowSpeeds = speeds(windowStart:windowEnd);
%
%         medianVelocity = median(windowSpeeds);
%         velocityThreshold = medianVelocity * params.velocityThresholdMultiplier;
%
%         if speeds(i) > velocityThreshold
%             saccadeCandidates(i) = 1;
%         end
%     end
%
%
%     % 3. Clustering
%     saccades = [];
%     inSaccade = false;
%     saccadeStart = 0;
%
%     for i = 1:length(saccadeCandidates)
%         if saccadeCandidates(i) && ~inSaccade % Saccade start
%             inSaccade = true;
%             saccadeStart = i;
%         elseif ~saccadeCandidates(i) && inSaccade % Saccade end
%             inSaccade = false;
%             saccadeEnd = i - 1;
%             saccades = [saccades; saccadeStart, saccadeEnd];
%         end
%     end
%
%
%         %Catch any saccade ending at the end of the data
%     if inSaccade
%         saccadeEnd = length(saccadeCandidates);
%         saccades = [saccades; saccadeStart, saccadeEnd];
%     end
%
%
%     % 4. Filtering
%     [saccades , RMSdist] = filterSaccades(eyePositions, saccades, samplingRate, params.minDuration, params.minAmplitude);
% end


function [saccades, RMSdist] = detectSaccadesORG(eyePositions, samplingRate, params)
% DETECTSACCADES Detects saccades in eye position data with a loss of eye
% metric.
%
%   Inputs:
%       eyePositions: 2D array (n x 2) of eye positions (x, y).
%       samplingRate: Sampling rate of the data (Hz).
%       params: Structure containing detection parameters:
%           .velocityThreshold:  Velocity threshold multiplier relative to median velocity
%           .minDuration:  Minimum saccade duration (ms)
%           .minAmplitude: Minimum saccade amplitude (pixels)
%           .lossWindow: Time window in ms for the loss metric calculation.
%           .lossThreshold: Threshold for the loss metric variance value.
%
%   Outputs:
%       saccades: A matrix where each row is a saccade, with columns:
%              start index, end index.
%       isLost: A boolean vector where a 1 indicates a point might be bad data due to loss of the eye.
% Parameter handling (set defaults if no params passed)
if nargin < 3 || isempty(params)
    params.velocityThreshold = 4;
    params.minDuration = 20; %ms
    params.minAmplitude = 0.5; %pixels
    params.lossWindow = 100; %ms
    params.lossThreshold = 2; %variance of 2 pixels^2
    warning('Parameter structure not found or empty; setting default saccade detection parameters.');
end


% 1. Calculate Velocities
velocities = calculateVelocity(eyePositions, samplingRate);
speeds = sqrt(sum(velocities.^2, 2)); % Magnitude of the velocities


% 2. Thresholding
medianVelocity = median(speeds);
velocityThreshold = medianVelocity * params.velocityThreshold; % dynamic threshold based on a multiple of the median
saccadeCandidates = speeds > velocityThreshold;

% 3. Clustering
saccades = [];
inSaccade = false;
saccadeStart = 0;

for i = 1:length(saccadeCandidates)
    if saccadeCandidates(i) && ~inSaccade % Saccade start
        inSaccade = true;
        saccadeStart = i;
    elseif ~saccadeCandidates(i) && inSaccade % Saccade end
        inSaccade = false;
        saccadeEnd = i - 1;
        saccades = [saccades; saccadeStart, saccadeEnd];
    end
end

%Catch any saccade ending at the end of the data
if inSaccade
    saccadeEnd = length(saccadeCandidates);
    saccades = [saccades; saccadeStart, saccadeEnd];
end

% 5. Filter using  loss metric and duration/amplitude
[saccades , RMSdist] = filterSaccades(eyePositions, saccades,...
    samplingRate, params.minDuration, params.minAmplitude);

% 6. Output data
% velocities
end



function velocities = calculateVelocity(positions, samplingRate)
% Calculates the velocities from position data.
% Assumes positions are in the format [x_pos, y_pos]
% Returns velocity at each timepoint (except for the first one).
% The final velocity will be the same as the penultimate velocity.
%
diffs = diff(positions);
velocities = diffs * samplingRate;
% Velocity at the last position will be the same as the velocity for
% the penultimate position
velocities = [velocities; velocities(end,:)];
end


function [filteredSaccades,RMSdist] = filterSaccades(eyePositions, saccades, samplingRate, minDuration, minAmplitude)
% Filters saccade candidates based on minimum duration, amplitude and data spread.
filteredSaccades = [];
RMSdist = [];
for i = 1:size(saccades, 1)
    startIndex = saccades(i, 1);
    endIndex = saccades(i, 2);

    % Basic error handling
    if startIndex < 1 || endIndex > size(eyePositions, 1) || startIndex >= endIndex
        fprintf('Saccade candidate %d has invalid start or end indices.\n', i);
        continue; % Skip to the next saccade
    end
    saccadePositions = eyePositions(startIndex:endIndex, :);

    % 1. Duration Check
    durationMs = (endIndex - startIndex + 1) / samplingRate * 1000;
    if durationMs < minDuration
        fprintf('Saccade candidate %d failed duration test. Duration = %.2f ms.\n', i, durationMs);
        continue; % Skip to the next check
    end

    % 2. Amplitude Check
    amplitude = sqrt(sum((saccadePositions(end,:) - saccadePositions(1,:)).^2)); % Euclidian distance between start and end
    if amplitude < minAmplitude
        fprintf('Saccade candidate %d failed amplitude test. Amplitude = %.2f pixels.\n', i, amplitude);
        continue; % Skip to the next check
    end

    % 3. Data Spread Check
    % If saccade contains data loss (any point is flagged), discard
    distances = [];
    %Calculate Euclidean distance between each pair of points
    numPoints = size(saccadePositions,1);
    for ii = 1:numPoints
        if ii == numPoints
            continue
        else
            dist = sqrt(sum((saccadePositions(ii,:) - saccadePositions(ii+1,:)).^2));
            distances = [distances; dist];
        end
    end

    %Compute the RMS value
    rmsDistance = sqrt(mean(distances.^2));

    if rmsDistance > 45
        fprintf('Saccade candidate %d failed rms test. RMS = %.2f pixels.\n', i, rmsDistance);
        continue
    end

    filteredSaccades = [filteredSaccades; saccades(i, :)];
    RMSdist = [RMSdist ; rmsDistance];

end
end



function saccadeMetrics = computeSaccadeMetrics(eyePositions, saccades, samplingRate)
%COMPUTESACCADEMETRICS Computes various metrics for detected saccades.
%
%   Inputs:
%       eyePositions: 2D array (n x 2) of eye positions (x, y).
%       saccades: A matrix where each row is a saccade, with columns:
%              start index, end index.
%       samplingRate: Sampling rate of the data (Hz).
%
%   Outputs:
%       saccadeMetrics: A struct where each field is a metric, as
%         follows:
%              meanVelocity:  Mean velocity for each saccade (pixels/sec).
%              peakVelocity:  Peak velocity for each saccade (pixels/sec).
%              direction:  Direction of each saccade (degrees).
%              amplitude:  Amplitude of each saccade (pixels).
%
numSaccades = size(saccades, 1);
saccadeMetrics.meanVelocity = zeros(numSaccades, 1);
saccadeMetrics.peakVelocity = zeros(numSaccades, 1);
saccadeMetrics.direction = zeros(numSaccades, 1);
saccadeMetrics.amplitude = zeros(numSaccades, 1);

for i = 1:numSaccades
    startIndex = saccades(i, 1);
    endIndex = saccades(i, 2);

    % Basic error handling
    if startIndex < 1 || endIndex > size(eyePositions, 1) || startIndex >= endIndex
        fprintf('Saccade %d has invalid start or end indices.\n', i);
        continue; % Skip to the next saccade
    end


    saccadePositions = eyePositions(startIndex:endIndex, :);

    % 1. Calculate Velocities
    % [avgVelocity, peakVelocity] = calculateVelocity2(saccadePositions);
    velocities = calculateVelocity(saccadePositions, samplingRate);
    speeds = sqrt(sum(velocities.^2, 2)); % Magnitude of the velocities

    % NEED to determine scaling factor -- default now to 10.
    speeds = round(speeds)/10;

    % 2. Mean Velocity
    saccadeMetrics.meanVelocity(i) = mean(speeds);

    % 3. Peak Velocity
    saccadeMetrics.peakVelocity(i) = max(speeds);

    % 4. Amplitude (Euclidian distance between start and end)
    amplitude = sqrt(sum((saccadePositions(end,:) - saccadePositions(1,:)).^2));
    saccadeMetrics.amplitude(i) = amplitude;

    % 5. Direction (using arctan2)
    deltaPos = saccadePositions(end,:) - saccadePositions(1,:); % x,y change in positions
    saccadeMetrics.direction(i) = rad2deg(atan2(deltaPos(2),deltaPos(1))); %Angle in degrees

    % 6. duration
    saccadeMetrics.duration(i) = height(saccadePositions);
end
end



% function [avgVelocity, peakVelocity] = calculateVelocity2(data)
% % Calculate velocities between consecutive points
% velocities = sqrt(diff(data(:,1)).^2 + diff(data(:,2)).^2) ./ diff(data(:,3));
% avgVelocity = mean(velocities);
% peakVelocity = max(velocities);
% end
