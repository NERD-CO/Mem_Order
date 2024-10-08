function [] = eyeTRACKproc_PupilSize_MO(cleanedDataLOC, saveLOC, ptID)

% Set data paths
if ~exist(cleanedDataLOC,'dir')
    mkdir(cleanedDataLOC)
end

% [For every file in eyeDATA, loop through to run analyses?]
% save cleaned files in new folder to only have raw data in this folder
% need to cd to new folder to save

% NO SHORTEN - because NOW Table INCLUDES time MARKERS

% CD to mainPath for raw data
cd(saveLOC);

% % Get contents of eyeDATA as a variable
% test = dir('*.mat');
% fileList = {test.name};

% for loop to loop throguh fileList
%for fileS = 1:length(fileList)
eyeData_pt = append('eyeData_', ptID,'.mat');
tempFile_name = eyeData_pt;

% Set tempFile_name
saveFile_name1 = ['cl_', tempFile_name];

% Load in file
load(tempFile_name, 'outInfo');

% Set tempEye- loop through every eye in every variant within variantS
% Go into variantS and determine # variants there
sessSnum = length(fieldnames(outInfo));

% Other option: Extract field names of variantS ahead of time
sessSfieldN = fieldnames(outInfo);
% Extract field names of each variant in variantS

for i = 1:sessSnum

    switch sessSfieldN{i}
        case 'encoding'
            % curTTLtable = outInfo.(sessSfieldN{i}).TTLinfo;
            % Clean up check
            % curTTLtable = curTTLtable(cellfun(@(x) ~isempty(x), curTTLtable , 'UniformOutput',true));
            outInfo.(sessSfieldN{i}).TTLinfo = outInfo.(sessSfieldN{i}).TTLinfo(1:90);
        case 'sceneRecog'
            outInfo.(sessSfieldN{i}).TTLinfo = outInfo.(sessSfieldN{i}).TTLinfo(1:180);
        case 'timeD'
            outInfo.(sessSfieldN{i}).TTLinfo = outInfo.(sessSfieldN{i}).TTLinfo(1:90);
    end

    % Loop through each eye in session
    for eyE = 1:2

        switch eyE
            case 1 % LEFT EYE
                % Get mean
                outInfo.(sessSfieldN{i}).leftEYE.oT_pupilS_meanCl =...
                    cleanPSmean(outInfo.(sessSfieldN{i}).leftEYE.oT_pupilS_mean);
                % Insert Nans in missing pupil values
                outInfo.(sessSfieldN{i}).leftEYE.oT_pupilS_rawCL =...
                    cleanEyeProcBF(outInfo.(sessSfieldN{i}).leftEYE.oT_pupilS_raw,0,...
                    outInfo.(sessSfieldN{i}).TTLinfo);

            case 2 % RIGHT EYE
                % Get mean
                outInfo.(sessSfieldN{i}).rightEYE.oT_pupilS_meanCl =...
                    cleanPSmean(outInfo.(sessSfieldN{i}).rightEYE.oT_pupilS_mean);

                % Insert Nans in missing pupil values
                outInfo.(sessSfieldN{i}).rightEYE.oT_pupilS_rawCL =...
                    cleanEyeProcBF(outInfo.(sessSfieldN{i}).rightEYE.oT_pupilS_raw,0,...
                    outInfo.(sessSfieldN{i}).TTLinfo);
        end % End of switch case
    end % End of Eye Loop
end % End of Session Loop
cd(cleanedDataLOC);
save(saveFile_name1, 'outInfo','-v7.3');


end



function [outMean] = cleanPSmean(inMean)

cleanPupil_size = inMean < 120; %create logical vector
inMean(cleanPupil_size) = nan;
outMean = inMean;

end


function [finalACSI] = cleanEyeProcBF(oldRAW,plotUi,TTLtable)


% PRE STEP 1 - SUPER outlier TRIALS not single events
oldRAW_SO = cell2mat(oldRAW);
oldRmean = mean(oldRAW_SO,'omitnan');
oldRstd = std(oldRAW_SO,'omitnan');
% oldRrms = rms(oldRAW_SO);
oldUthreh = oldRmean + (oldRstd * 2);
oldDthreh = 25;

% STACK EVERY TRIAL within EYE / CONDITION

for tsI = 1:height(oldRAW)

    tmpTrace = oldRAW{tsI};
    tmpIndex = tmpTrace < oldDthreh | tmpTrace > oldUthreh;
    tmpTrace(tmpIndex) = nan;

    oldRAW{tsI} = tmpTrace;

end


% PRE STEP 2 - SUPER outlier TRIALS IN LENGTH
for tsI = 1:height(oldRAW)

    disp(['Trial analysis begin ', num2str(tsI)])

    tmpTrace = oldRAW{tsI};
    if isempty(tmpTrace)
        oldRAW{tsI} = NaN;
        continue
    else
        tmpTTL = TTLtable{tsI}.ELNKint(5);

        if tmpTTL > 9000
            tmpTrace(5001:end) = nan;
        else
            tmpTrace(tmpTTL + 500:end) = nan;
        end
        oldRAW{tsI} = tmpTrace;
    end
end


% figure;
% for nfi = 1:length(oldRAW)
%     tmpEFi = oldRAW{nfi};
%     plot(tmpEFi,'Color',[0.3 0.3 0.3 0.3],'LineWidth',1.5)
%     hold on
%     xline(TTLtable{nfi}.ELNKint(5))
% end
% title('Remove OUTlier trials')
% % pause(10)
% close all



% STEP 1
trialIDS = zeros(height(oldRAW),1,'logical');
numCpoints = zeros(height(oldRAW),1);
locCpoints = cell(height(oldRAW),1);
for tsI = 1:height(oldRAW)

    tmpTrace = oldRAW{tsI};
    tmpTUpth = mean(tmpTrace, 'omitnan') + (std(tmpTrace, 'omitnan')*2);
    tmpTDpth = mean(tmpTrace, 'omitnan') - (std(tmpTrace, 'omitnan')*2);

    tmpIndex = tmpTrace < tmpTDpth | tmpTrace > tmpTUpth; % FIGURE OUT UPPER LIMIT??????????

    % plot(tmpTrace); yline(tmpTUpth); yline(tmpTDpth)
    % pause 
    % cla

    if sum(tmpIndex) == 0
        continue
    else
        changeIndices = find(diff(tmpIndex));
        numCpoints(tsI) = sum(diff(tmpIndex) ~= 0);
        locCpoints{tsI} = changeIndices;
        trialIDS(tsI) = true;
        % FOR PLOT CHECKING................................................
        if plotUi
            plot(tmpTrace);
            hold on
            plot(tmpIndex*100)
            plot(changeIndices,ones(numel(changeIndices),1)*500,'k.','MarkerSize',25)
            xline(locCpoints{tsI})
            pause
            cla
        end
        % FOR PLOT CHECKING................................................
    end
end

% STEP 2
trialIDSf = find(trialIDS);
numCptsf = numCpoints(trialIDS);
locCptsf = locCpoints(trialIDS);

% Loop through each affected trial
forwardPts = cell(height(trialIDSf),1);
backwardPts = cell(height(trialIDSf),1);
for efi = 1:length(trialIDSf)

    tmpEFi = oldRAW{trialIDSf(efi)};

    % Loop through change points
    for chpi = 1:numCptsf(efi)

        changLoc = locCptsf{efi}(chpi);

        % forward point of change
        tmpWindowStart = changLoc - 5;
        tmpWindowStop = changLoc;
        if tmpWindowStart < 1
            forwardPOINT = 1; % consider tmpNRMS block of code to find index
        else

            chek1 = true;
            while chek1

                if tmpWindowStart < 1
                    forwardPOINT = 1;
                    chek1 = false;
                else
                    tmpVecC = tmpEFi(tmpWindowStart:tmpWindowStop);
                    tmpNRMS = abs(tmpVecC/rms(tmpEFi(tmpEFi > 100))-1);

                    if ~any(tmpNRMS < 0.1) % all changed
                        tmpWindowStop = tmpWindowStart;
                        tmpWindowStart = tmpWindowStart - 5;
                        continue
                    else
                        rmsChangpt = find(tmpNRMS < 0.1,1,'first');
                        efIindexVec = tmpWindowStart:tmpWindowStop;
                        forwardPOINT = efIindexVec(rmsChangpt);
                        chek1 = false;
                    end
                end
            end
        end % END of FORWARD IF/ELSE

        % backward point of change
        tmpWindowStart = changLoc;
        tmpWindowStop = changLoc + 5;
        if tmpWindowStop > numel(tmpEFi)
            backwardPOINT = changLoc;
        else
            chek2 = true;
            while chek2

                if tmpWindowStop > numel(tmpEFi)
                    backwardPOINT = numel(tmpEFi);
                    chek2 = false;
                else
                    tmpVecC = tmpEFi(tmpWindowStart:tmpWindowStop);
                    tmpNRMS = abs(tmpVecC/rms(tmpEFi(tmpEFi > 100))-1);

                    if ~any(tmpNRMS < 0.1) % all changed
                        tmpWindowStart = tmpWindowStop;
                        tmpWindowStop = tmpWindowStart + 5;
                        continue
                    else
                        rmsChangpt = find(tmpNRMS < 0.1,1,'first');
                        efIindexVec = tmpWindowStart:tmpWindowStop;
                        backwardPOINT = efIindexVec(rmsChangpt);
                        chek2 = false;
                    end
                end

            end % END of WHILE LOOP for BACKWARD
        end % END of BACKWARD IF/ELSE
        forwardPts{efi}(chpi) = forwardPOINT;
        backwardPts{efi}(chpi) = backwardPOINT;
    end % END of CHANGE POINT LOop

end

% STEP 3 - INSERT NANS
newRAW = oldRAW;

for naNdi = 1:height(trialIDSf)
    newVecEfi = newRAW{trialIDSf(naNdi)};
    for inNan = 1:numel(forwardPts{naNdi})
        nanSTART = forwardPts{naNdi}(inNan);
        nanSTOP = backwardPts{naNdi}(inNan);
        newVecEfi(nanSTART:nanSTOP) = nan;
    end
    newRAW{trialIDSf(naNdi)} = newVecEfi;
end

% PLOT STEP 3 
% figure;
% for nfi = 1:length(newRAW)
%     tmpEFi = newRAW{nfi};
%     plot(tmpEFi,'Color',[0.3 0.3 0.3 0.3],'LineWidth',1.5)
%     hold on
% end
% title('JAT Clean up')
% close all

% STEP 4 FIRST DERIVATIVE MASK
dernewRAW = newRAW;
for drIi = 1:length(newRAW)

    tmpNEW = newRAW{drIi};
    % First derivative of x
    dx = diff(tmpNEW);
    % Standard deviation of dx
    sigma = std(dx,'omitnan');
    % 2.5 times the standard deviation
    result = 2.5 * sigma;
    absDer = abs(dx);
    firDerSDmask = absDer > result;

    derNEW = tmpNEW;
    derNEW(firDerSDmask) = nan;

    dernewRAW{drIi} = derNEW;
end

% PLOT STEP 4 
% figure;
% for nfi = 1:length(dernewRAW)
%     tmpEFi = dernewRAW{nfi};
%     plot(tmpEFi,'Color',[0.3 0.3 0.3 0.3],'LineWidth',1.5)
%     hold on
% end
% title('1st Deriv')
% close all

% STEP 5
% 3rd ORDER ButterWORTH with 4 Hz cut off
% Define the filter order
order = 3;
% Define the cutoff frequency
fc = 4; % Hz
% Define the sampling frequency
fs = 1000; % Hz
% Normalize the cutoff frequency to the Nyquist frequency
Wn = fc / (fs/2);
% Design the filter coefficients
[b, a] = butter(order, Wn);
% Filter the signal using the filter coefficients

fourHzbutter = dernewRAW;
for buIi = 1:length(dernewRAW)
    tmpNEWt = dernewRAW{buIi};
    % Interpolate NaNs in the data
    nanLOCS = isnan(tmpNEWt);
    t = 1:length(tmpNEWt);

    % FOR ERROR CHECKING
    % disp(buIi)

    allnanCheck = sum(isnan(tmpNEWt))/numel(tmpNEWt);
    if allnanCheck > 0.85
        fourHzbutter{buIi} = transpose(tmpNEWt);
        continue
    end

    data_interp1 = interp1(t(~isnan(tmpNEWt)), tmpNEWt(~isnan(tmpNEWt)), t);
    if any(isnan(data_interp1))
        dataMEAN = mean(data_interp1,'omitnan');
        data_interp1(isnan(data_interp1)) = dataMEAN;
    end
    % Apply the filter to the interpolated data
    filtered_data = filtfilt(b, a, data_interp1);
    % filtered_signal = filtfilt(b, a, tmpNEWt);
    fbTfil = filtered_data;
    % PUT NANS back
    fbTfil(nanLOCS) = nan;
    fourHzbutter{buIi} = fbTfil;
end

finalNEWot = fourHzbutter;


% PLOT STEP 5 
figure;
for nfi = 1:length(finalNEWot)
    tmpEFi = finalNEWot{nfi};
    plot(tmpEFi,'Color',[0.3 0.3 0.3 0.3],'LineWidth',1.5)
    hold on
end
title('4Hz low cut')
close all
maxLEN = max(cellfun(@(x) numel(x), finalNEWot, 'UniformOutput',true));
nanMAT = nan(length(finalNEWot),maxLEN);
for ni = 1:length(finalNEWot)
    curLEN = numel(finalNEWot{ni});
    nanMAT(ni,1:curLEN) = finalNEWot{ni};
end


% FIND last NAN
firstAllNAN = find(sum(isnan(nanMAT),1) == length(finalNEWot),1,'first');

nanMATn = nanMAT(:,1:firstAllNAN);

fracNAN = zeros(1,width(nanMATn));
for ini = 1:width(nanMATn)

    fracNAN(ini) = sum(isnan(nanMATn(:,ini)))/height(nanMATn);

end

columnTrim = find(fracNAN >= .80, 1, 'first');

nanMATn2 = nanMATn(:,1:columnTrim);


[corVals] = computeASCI_eye_v2(nanMATn2);
% ccCUToff = quantile(corVals,0.12);
if std(corVals,'omitnan') < 0.29
    ccCUToff = mean(corVals,'omitnan') - (std(corVals,'omitnan')*1.65);
else
    ccCUToff = mean(corVals,'omitnan') - (std(corVals,'omitnan')*0.25);
end
corInds = corVals > ccCUToff;
% corInds = corVals > 0.6;

finalACSI = finalNEWot;
finalACSI(~corInds) = {nan};


% PLOT STEP 4 
% nanMATntemp = nanMATn(corInds,:);
% figure;
% for nfi = 1:height(nanMATntemp)
%     tmpEFi = nanMATntemp(nfi,:);
%     plot(tmpEFi,'Color',[0.3 0.3 0.3 0.3],'LineWidth',1.5)
%     hold on
% end
% title('1st Deriv')
% close all


% SMOOTH GUASSION or SANGOLY

% finalSMOOTH = fourHzbutter;
% for smIi = 1:length(fourHzbutter)
% 
%     tmpNEWts = fourHzbutter{smIi};
% 
%     finSMdata = smoothdata(tmpNEWts , "sgolay" , 60);
% 
%     finalSMOOTH{smIi} = finSMdata;
% end
% 
% figure;
% for nfi = 1:length(finalSMOOTH)
%     tmpEFi = finalSMOOTH{nfi};
%     plot(tmpEFi,'Color',[0.3 0.3 0.3 0.3],'LineWidth',1.5)
%     hold on
% end
% title('Final smooth')



end




% function [newTRIMvec] = shortenVEC(oldTRIMvec)
% 
% % Find shortest row in column 6
% min_varEye = min(cellfun(@(x) numel(x), oldTRIMvec, 'UniformOutput', true));
% % Trim all rows in col 6 to shortest length
% newTRIMvec = oldTRIMvec;
% for trI = 1:height(oldTRIMvec)
%     newTRIMvec{trI} = oldTRIMvec{trI}(1:min_varEye);
% end
% 
% 
% end

