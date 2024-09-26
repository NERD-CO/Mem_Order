function [xBounds , yBounds] = getPICtureDims_MO(imageLOCATION,imageNAME)

[~,~,tmpEXT] = fileparts(imageNAME);

switch tmpEXT
    case '.mp4'
        cd(imageLOCATION)
        videoObj = VideoReader([imageLOCATION filesep imageNAME]);

        % Get the image dimensions
        imageWidth = videoObj.Width;
        imageHeight = videoObj.Height;

    case '.jpeg'

        img = imread(imageLOCATION);
        [xRow,yCol] = size(img,[1,2]);

        imageWidth = xRow;
        imageHeight = yCol;

end

screenWidth = 1920;
screenHeight = 1080;

% Calculate centered coordinates -- top left corner of centered image
x = (screenWidth - imageWidth) / 2;
y = (screenHeight - imageHeight) / 2;

xBounds = round([x x+imageWidth]);
yBounds = round([y y+imageHeight]);

end