function [xBounds , yBounds] = getPICtureDims_MO(imageLOCATION,imageNAME)

[~,~,tmpEXT] = fileparts(imageNAME);

switch tmpEXT
    case '.mp4'
        cd(imageLOCATION)
        videoObj = VideoReader(imageNAME);

        % Get the image dimensions
        frameWidth = videoObj.Width;
        frameHeight = videoObj.Height;


    case '.jpeg'


end

img = imread(imageLOCATION);
[xRow,yCol] = size(img,[1,2]);
screenWidth = 1024;
screenHeight = 768;
% Image dimensions
imageWidth = xRow;
imageHeight = yCol;

% Calculate centered coordinates -- top left corner of centered image
x = (screenWidth - imageWidth) / 2;
y = (screenHeight - imageHeight) / 2;

xBounds = round([x x+imageWidth]);
yBounds = round([y y+imageHeight]);

end