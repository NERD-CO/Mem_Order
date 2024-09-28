function [outEye1 , outEye2] = createEYEtable_MO_GAZE(pos1,pos2,posF1,posF2,...
    posP1,posP2,timeIN)

for ei = 1:2
    switch ei
        case 1
            eyedataR = {pos1};
            eyeCenR = {mean(pos1,'omitnan')};
            eyeSDR = {std(pos1,'omitnan')};

            eyedataF = {pos1(posF1,:)};
            eyeCenF = {mean(pos1(posF1,:),'omitnan')};
            eyeSDF = {std(pos1(posF1,:),'omitnan')};

            eyedataP = {pos1(posP1,:)};
            eyeCenP = {mean(pos1(posP1,:),'omitnan')};
            eyeSDP = {std(pos1(posP1,:),'omitnan')};

            timeFr = {timeIN(posF1)};
            timePic = {timeIN(posP1)};

            outEye1 = table(eyedataR,eyeCenR,eyeSDR,...
                eyedataF,eyeCenF,eyeSDF,...
                eyedataP,eyeCenP,eyeSDP,...
                timeFr,timePic,'VariableNames',...
                {'gaze_Raw','gaze_Raw_cen','gaze_Raw_sd',...
                'gaze_Scr','gaze_Scr_cen','gaze_Scr_sd',...
                'gaze_Pic','gaze_Pic_cen','gaze_Pic_sd',...
                'gaze_Scr_time','gaze_Pic_time'});

        case 2

            eyedataR = {pos2};
            eyeCenR = {mean(pos2,'omitnan')};
            eyeSDR = {std(pos2,'omitnan')};

            eyedataF = {pos2(posF2,:)};
            eyeCenF = {mean(pos2(posF2,:),'omitnan')};
            eyeSDF = {std(pos2(posF2,:),'omitnan')};

            eyedataP = {pos2(posP2,:)};
            eyeCenP = {mean(pos2(posP2,:),'omitnan')};
            eyeSDP = {std(pos2(posP2,:),'omitnan')};

            timeFr = {timeIN(posF2)};
            timePic = {timeIN(posP2)};
 
            outEye2 = table(eyedataR,eyeCenR,eyeSDR,...
                eyedataF,eyeCenF,eyeSDF,...
                eyedataP,eyeCenP,eyeSDP,...
                timeFr,timePic,'VariableNames',...
                {'gaze_Raw','gaze_Raw_cen','gaze_Raw_sd',...
                'gaze_Scr','gaze_Scr_cen','gaze_Scr_sd',...
                'gaze_Pic','gaze_Pic_cen','gaze_Pic_sd',...
                'gaze_Scr_time','gaze_Pic_time'});

    end


end
end
