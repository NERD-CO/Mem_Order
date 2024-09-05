function [outEye1 , outEye2] = createEYEtable(ps1,ps2,pos1,pos2)

%% function [outEye1 , outEye2] = createEYEtable(ps1,ps2,pos1,pos2)

%% Loop over eyes
for ei = 1:2
    switch ei
        case 1
            eyedata = {pos1};
            eyeCen = {mean(pos1)};
            eyeSD = {std(pos1)};
            Q3pos = quantile(pos1,0.75);
            Q1pos = quantile(pos1,0.25);
            eyeCD = (Q3pos - Q1pos) / (Q3pos + Q1pos);
            eyedist = {pdist2(eyeCen{1},eyedata{1},'euclidean')};
            pupdata = {ps1};
            pupCen = mean(ps1);
            pupSD = std(ps1);
            Q3pup = quantile(ps1,0.75);
            Q1pup = quantile(ps1,0.25);
            pupCD = (Q3pup - Q1pup) / (Q3pup + Q1pup);

            outEye1 = table(eyedata,eyeCen,eyeSD,eyeCD,eyedist,...
                pupdata,pupCen,pupSD,pupCD,'VariableNames',...
                {'oT_posit_raw','oT_posit_cen','oT_posit_sd','oT_posit_cd',...
                'oT_posit_dist','oT_pupilS_raw','oT_pupilS_mean','oT_pupilS_sd',...
                'oT_pupilS_cd'});
        case 2
            eyedata = {pos2};
            eyeCen = {mean(pos2)};
            eyeSD = {std(pos2)};
            Q3pos = quantile(pos2,0.75);
            Q1pos = quantile(pos2,0.25);
            eyeCD = (Q3pos - Q1pos) / (Q3pos + Q1pos);
            eyedist = {pdist2(eyeCen{1},eyedata{1},'euclidean')};
            pupdata = {ps2};
            pupCen = mean(ps2);
            pupSD = std(ps2);
            Q3pup = quantile(ps2,0.75);
            Q1pup = quantile(ps2,0.25);
            pupCD = (Q3pup - Q1pup) / (Q3pup + Q1pup);

            outEye2 = table(eyedata,eyeCen,eyeSD,eyeCD,eyedist,...
                pupdata,pupCen,pupSD,pupCD,'VariableNames',...
                {'oT_posit_raw','oT_posit_cen','oT_posit_sd','oT_posit_cd',...
                'oT_posit_dist','oT_pupilS_raw','oT_pupilS_mean','oT_pupilS_sd',...
                'oT_pupilS_cd'});
    end
end

end