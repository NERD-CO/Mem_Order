function [outfiles] = getfiles(dirIN,stage,ftype)

    %% Go to directory
    cd(dirIN);

    %% Get files
    switch stage
        case 1
            foldeS = dir();
            foldeS2 = {foldeS.name};
            foldeS3 = foldeS2(~ismember(foldeS2,{'.','..'}));
            outfiles = foldeS3;

        case 2
            filES = dir(['*.',ftype]);
            filES2 = {filES.name};
            outfiles = filES2;
    end
end