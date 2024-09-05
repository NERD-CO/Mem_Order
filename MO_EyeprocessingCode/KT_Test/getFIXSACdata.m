function [fixTABLE_Eye0,fixTABLE_Eye1,sacTABLE_Eye0,sacTABLE_Eye1] = getFIXSACdata(eyesOUT,gxALL,gyALL,timeALL)

%% function [fixTABLE_Eye0,fixTABLE_Eye1,sacTABLE_Eye0,sacTABLE_Eye1] = getFIXSACdata(eyesOUT,gxALL,gyALL,timeALL)

%% Fixation extract
eyeIDs = [0 1];
fix_eyeZEROp = cell(height(eyesOUT.eye_0.Fixes),1);
fix_eyeZEROd = nan(height(eyesOUT.eye_0.Fixes),1);
fix_eyeONEp = cell(height(eyesOUT.eye_1.Fixes),1);
fix_eyeONEd = nan(height(eyesOUT.eye_1.Fixes),1);

%% Loop over eyes
for eye2u = 1:2
    tmpEYE = eyeIDs(eye2u);
    tmpFIXES = eyesOUT.(['eye_',num2str(tmpEYE)]).Fixes;

    %% Loop over fixes
    for fixI = 1:height(tmpFIXES)
        startINDf = tmpFIXES.StartTime(fixI);
        endINDf = tmpFIXES.EndTime(fixI);

        startTTf = find(timeALL == startINDf);
        endTTf = find(timeALL == endINDf);

        switch tmpEYE
            case 0
                fix_eyeZEROp{fixI} = array2table([transpose(gxALL(1,startTTf:endTTf)) ...
                    transpose(gyALL(1,startTTf:endTTf))...
                    transpose(timeALL(1,startTTf:endTTf))],...
                    'VariableNames',{'gx','gy','timeTTL'});
                fix_eyeZEROd(fixI) = (endTTf - startTTf)/1000;
            case 1
                fix_eyeONEp{fixI} =...
                    array2table([transpose(gxALL(1,startTTf:endTTf))...
                    transpose(gyALL(1,startTTf:endTTf))...
                    transpose(timeALL(1,startTTf:endTTf))],...
                    'VariableNames',{'gx','gy','timeTTL'});
                fix_eyeONEd(fixI) = (endTTf - startTTf)/1000;
        end
    end
end

%% Fix eye tables
fixTABLE_Eye0 = table(fix_eyeZEROp , fix_eyeZEROd,...
    'VariableNames',{'FIX_Eye0_Points','FIX_Eye0_DurSecs'});

fixTABLE_Eye1 = table(fix_eyeONEp , fix_eyeONEd,...
    'VariableNames',{'FIX_Eye1_Points','FIX_Eye1_DurSecs'});

%% Preallocate
sac_eyeZEROp = cell(height(eyesOUT.eye_0.Saccades),1);
sac_eyeZEROd = nan(height(eyesOUT.eye_0.Saccades),1);
sac_eyeONEp = cell(height(eyesOUT.eye_1.Saccades),1);
sac_eyeONEd = nan(height(eyesOUT.eye_1.Saccades),1);

%% Loop over eyes
for eye2u = 1:2
    tmpEYE = eyeIDs(eye2u);
    tmpSACCADES = eyesOUT.(['eye_',num2str(tmpEYE)]).Saccades;

    %% Loop over saccades
    for sacI = 1:height(tmpSACCADES)
        startINDs = tmpSACCADES.StartTime(sacI);
        endINDs = tmpSACCADES.EndTime(sacI);

        startTTs = find(timeALL == startINDs);
        endTTs = find(timeALL == endINDs);

        switch tmpEYE
            case 0
                sac_eyeZEROp{sacI} =...
                    array2table([transpose(gxALL(1,startTTs:endTTs))...
                    transpose(gyALL(1,startTTs:endTTs))...
                    transpose(timeALL(1,startTTs:endTTs))],...
                    'VariableNames',{'gx','gy','timeTTL'});
                sac_eyeZEROd(sacI) = (endTTs - startTTs)/1000;
            case 1
                sac_eyeONEp{sacI} =...
                    array2table([transpose(gxALL(1,startTTs:endTTs))...
                    transpose(gyALL(1,startTTs:endTTs))...
                    transpose(timeALL(1,startTTs:endTTs))],...
                    'VariableNames',{'gx','gy','timeTTL'});
                sac_eyeONEd(sacI) = (endTTs - startTTs)/1000;
        end
    end
end

%% Fix saccades
sacTABLE_Eye0 = table(sac_eyeZEROp , sac_eyeZEROd,...
    'VariableNames',{'SAC_Eye0_Points','SAC_Eye0_DurSecs'});

sacTABLE_Eye1 = table(sac_eyeONEp , sac_eyeONEd,...
    'VariableNames',{'SAC_Eye1_Points','SAC_Eye1_DurSecs'});

end