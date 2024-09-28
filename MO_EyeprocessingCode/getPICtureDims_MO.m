function [xBounds , yBounds] = getPICtureDims_MO(imageLOCATION,...
    picTABLE , trialNUMBER , SessionTYPE)
% picTABLE, tttrialir, SessionTYPE

trialROW = picTABLE(trialNUMBER,:);

% [~,~,tmpEXT] = fileparts(imageNAME);

switch SessionTYPE
    case 'en'
        cd(imageLOCATION)

        imageNAME = trialROW.ClipName{1};  
        videoObj = VideoReader(imageNAME);
        % Get the image dimensions
        imageWidth = videoObj.Width;
        imageHeight = videoObj.Height;

    case 'sc'

        cd(imageLOCATION)

        imageNAME = trialROW.ClipName{1};
        img = imread(imageNAME);
        [xRow,yCol] = size(img,[1,2]);

        imageWidth = xRow;
        imageHeight = yCol;

    case 'ti'
        cd(imageLOCATION)

        % Get name from clipid
        imageNAME = trialROW.ClipName{1};
        % Get order from frame order
        imageORDER = trialROW.FrameORDER{1};
        % Create loader and stick together 
        fNAMEt = extractBefore(imageNAME,'.');
        exTEN = extractAfter(imageNAME,'.');

        imageONE = [fNAMEt(1:end-1),num2str(imageORDER(1)),'.',exTEN];
        imageTWO = [fNAMEt(1:end-1),num2str(imageORDER(2)),'.',exTEN];

        % load in left
        img1 = imread(imageONE);
        % get height and width
        [xRow1,yCol1] = size(img1,[1,2]);

        imageWidth1 = xRow1;
        imageHeight1 = yCol1;

        % load in right
        img2 = imread(imageTWO);
        % get height and width
        [xRow2,~] = size(img2,[1,2]);

        imageWidth2 = xRow2;
        % imageHeight2 = yCol2;

        % height will equal 1
        imageHeight  = imageHeight1*0.6;
        % width will equal +
        imageWidth = imageWidth1 + imageWidth2;

end

screenWidth = 1920;
screenHeight = 1180;

% Calculate centered coordinates -- top left corner of centered image
x = (screenWidth - imageWidth) / 2;
y = (screenHeight - imageHeight) / 2;

xBounds = round([x x+imageWidth]);
yBounds = round([y y+imageHeight]);

end