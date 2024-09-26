function [] = eyeTRACKproc_PupilLocation_MO(cleanedDataLOC, saveLOC, ptID)

% Set data paths
if ~exist(cleanedDataLOC,'dir')
    mkdir(cleanedDataLOC)
end

% NO SHORTEN - because NOW Table INCLUDES time MARKERS

% CD to mainPath for raw data
cd(saveLOC);

eyeData_pt = append('eyeData_', ptID,'.mat');
tempFile_name = eyeData_pt;

% Set tempFile_name
saveFile_name1 = ['cl_', tempFile_name];

% Load in file
load(tempFile_name, 'outInfo');

sessSnum = length(fieldnames(outInfo));

% Other option: Extract field names of variantS ahead of time
sessSfieldN = fieldnames(outInfo);

for i = 1:sessSnum
    % Loop through each eye in variant
    for eyE = 1:2
        disp(eyE)
        % FOR ERROR CHECKING
        % disp(eyE)

        switch eyE
            case 1
                gazeData = [variantS.(varSfieldN{i}).dataTable.Left_L_gaze_Pic,...
                    variantS.(varSfieldN{i}).dataTable.Left_L_gaze_Pic_time];

                curVariant.Left_L_GAZEcl = cleanEye_GAZE(gazeData,0);
            case 2
                gazeData = [variantS.(varSfieldN{i}).dataTable.Left_R_gaze_Pic,...
                    variantS.(varSfieldN{i}).dataTable.Left_R_gaze_Pic_time];

                curVariant.Left_R_GAZEcl = cleanEye_GAZE(gazeData,0);

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Save cleaned data into variantS
        variantS.(varSfieldN{i}).dataTable = curVariant;

    end

end
cd(cleanedDataLOC);
save(saveFile_name1, 'variantS');
%end



end

function [finalACSI] = cleanEye_GAZE(oldRAW,plotUi)


maxDispersion = 10;
minDuration = 75;

oldRAW1 = oldRAW;
finalACSI = cell(height(oldRAW1),1);
for oi = 1:height(oldRAW1)

    gazeDATA = [oldRAW1{oi,1},...
    double(oldRAW1{oi,2})];
    % Fixations
    fixations = detectFixations(gazeDATA, maxDispersion, minDuration);

    if plotUi
        [cenS , ~] = plotFixations(fixations , gazeDATA);
        xlim([0 1028])
        ylim([0 768])
    else
        [cenS , ~] = plotFixations(fixations , gazeDATA);
        close all
    end

    fixations1 = [fixations , cenS];

    if isempty(fixations)
        trialFixSacc.fixations = nan;
    else
        fixations2 = array2table(fixations1,'VariableNames',{'starttime','endtime',...
            'durationMS','fixCenX','fixCenY'});
        trialFixSacc.fixations = fixations2;
    end

    % detect saccades
    disp(oi)
    saccades = detectSaccades(gazeDATA, plotUi);
    if isnan(saccades)
        saccadeParams = nan;
    else
        saccadeParams = extractSaccadeParameters(gazeDATA, saccades);
        saccadeParams = struct2table(saccadeParams);
        saccadeParams.starttime = saccades(:,1);
        saccadeParams.endtime = saccades(:,2);
    end
    trialFixSacc.saccades = saccadeParams;
    finalACSI{oi} = trialFixSacc;

end

end





function fixations = detectFixations(gazeData, maxDispersion, minDuration)
    % Initialize variables
    fixations = [];
    startIndex = 1;
    numPoints = size(gazeData, 1);

    % Loop through gaze data points
    while startIndex < numPoints
        endIndex = startIndex;
        while endIndex <= numPoints && calculateDispersion(gazeData(startIndex:endIndex, :)) < maxDispersion
            endIndex = endIndex + 1;
        end
        endIndex = endIndex - 1;

        % Calculate duration
        duration = gazeData(endIndex, 3) - gazeData(startIndex, 3); % Assuming column 3 has timestamps
        if duration >= minDuration
            fixations = [fixations; gazeData(startIndex, 3), gazeData(endIndex, 3), duration];
        end

        % Move to the next possible start
        startIndex = endIndex + 1;
    end
end

function dispersion = calculateDispersion(points)
    % Calculate the maximum dispersion in x and y coordinates
    xCoords = points(:, 1); % Assuming column 1 has x coordinates
    yCoords = points(:, 2); % Assuming column 2 has y coordinates
    dispersion = max([max(xCoords) - min(xCoords), max(yCoords) - min(yCoords)]);
end



function [centroidsS , newCOLORs] = plotFixations(fixEpochs , gazeDATA)

newCOLORs = rand(height(fixEpochs),3);
centroidsS = zeros(height(fixEpochs),2);
for fii = 1:height(fixEpochs)

    fixTemp = fixEpochs(fii,:);

    fixStart = find(gazeDATA(:,3) == fixTemp(1));
    fixEnd = find(gazeDATA(:,3) == fixTemp(2));
    fixInd = fixStart:fixEnd;

    fixData = gazeDATA(fixInd,1:2);

    centroidsS(fii,:) = mean(fixData);

    hold on
    plot(fixData(:,1),fixData(:,2),'Color',newCOLORs(fii,:))

end


end




function saccades = detectSaccades(gazeData , plotUi)
    % Initialize variables
    gazeData = double(gazeData);
    numPoints = size(gazeData, 1);
    velocities = zeros(numPoints-1, 1);
    
    % Calculate velocities between consecutive points
    for i = 1:numPoints-1
        deltaX = gazeData(i+1, 1) - gazeData(i, 1);
        deltaY = gazeData(i+1, 2) - gazeData(i, 2);
        deltaT = gazeData(i+1, 3) - gazeData(i, 3);
        
        if deltaT == 0
            continue;  % Avoid division by zero
        end
        
        velocities(i) = sqrt(deltaX^2 + deltaY^2) / deltaT;
    end

    smVel = smoothdata(velocities,'gaussian',30);
    upvelThreshd = mean(smVel) + (std(smVel)*0.15);
    dnvelThreshd1 = mean(smVel) - (std(smVel)*0.15);

    % find above lower threshold
    close all
    smVel2 = smVel;
    smVel2(smVel < dnvelThreshd1) = mean(smVel(smVel < dnvelThreshd1));
    smVel2(1:25) = mean(smVel(smVel < dnvelThreshd1));
    % plot(smVel2);

    if length(velocities) < 40

        saccades = nan;
        return

    end

    [~,locs,~,~] = findpeaks(smVel2,'MinPeakHeight',upvelThreshd,'MinPeakDistance',30,...
        'Annotate','extents');

    if numel(locs) < 2

        saccades = nan;
        return

    end

    if plotUi

        findpeaks(smVel2,'MinPeakHeight',upvelThreshd,'MinPeakDistance',30,...
            'Annotate','extents');
  
    end

    startPeakInds = zeros(numel(locs),1);
    endPeakInds = zeros(numel(locs),1) ;
    for peakI = 1:length(locs)

        peakLOCi = locs(peakI);

        % dnvelThreshd1
        sii = peakLOCi;
        while smVel2(sii) > dnvelThreshd1
            sii = sii - 1;
        end
        startPeakInds(peakI) = sii;

        eii = peakLOCi;
        while smVel2(eii) > dnvelThreshd1
            eii = eii + 1;

            if eii > length(smVel2)
                eii = length(smVel2);
                break
            end

        end
        endPeakInds(peakI) = eii;

    end

    if plotUi

        xline(startPeakInds,'g-')
        xline(endPeakInds,'r-')
        legend('off')

        pause
        close all
    end

    % Detect saccades based on velocity threshold
    saccades = zeros(numel(startPeakInds),1);
    for i = 1:length(startPeakInds)
        saccades(i,1) = gazeData(startPeakInds(i), 3);
        saccades(i,2) = gazeData(endPeakInds(i), 3);
    end

end


function saccadeParams = extractSaccadeParameters(gazeData, saccades)
    % Input:
    % gazeData - Nx3 matrix [xPos, yPos, timestamps]
    % saccades - Mx2 matrix [startTimes, endTimes]
    
    % Output:
    % saccadeParams - structure array containing parameters for each saccade
    gazeData = double(gazeData);
    numSaccades = size(saccades, 1);
    saccadeParams(numSaccades) = struct();

    for i = 1:numSaccades
        startTime = saccades(i, 1);
        endTime = saccades(i, 2);
        
        % Get indices for the saccade
        startIndex = find(gazeData(:,3) == startTime, 1, 'first');
        endIndex = find(gazeData(:,3) == endTime, 1, 'last');
        
        % Extract saccade data
        saccadeData = gazeData(startIndex:endIndex, :);
        
        % Calculate parameters
        saccadeParams(i).amplitude = calculateAmplitude(saccadeData);
        saccadeParams(i).duration = saccadeData(end, 3) - saccadeData(1, 3);
        [saccadeParams(i).velocity, saccadeParams(i).peakVelocity] = calculateVelocity(saccadeData);
        saccadeParams(i).direction = calculateDirection(saccadeData);
    end
end

function amplitude = calculateAmplitude(data)
    % Calculate the Euclidean distance between the first and last points
    deltaX = data(end,1) - data(1,1);
    deltaY = data(end,2) - data(1,2);
    amplitude = sqrt(deltaX^2 + deltaY^2);
end

function [avgVelocity, peakVelocity] = calculateVelocity(data)
    % Calculate velocities between consecutive points
    velocities = sqrt(diff(data(:,1)).^2 + diff(data(:,2)).^2) ./ diff(data(:,3));
    avgVelocity = mean(velocities);
    peakVelocity = max(velocities);
end

function direction = calculateDirection(data)
    % Calculate the direction of the saccade in degrees
    deltaY = data(end,2) - data(1,2);
    deltaX = data(end,1) - data(1,1);
    direction = atan2d(deltaY, deltaX);
end

