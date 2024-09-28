heightScaler = 0.8; % 80% height of the screen y axis
imgHeight = screenYpixels * heightScaler;
imgWidth = imgHeight * orig_height / orig_width;
theRect = [0 0 imgWidth imgHeight];
dstRect = CenterRectOnPointd(theRect, screenXpixels/2, screenYpixels/2);

%%


% Get the size of the image
[orig_width, orig_height, ~]=size(img1);
% Define the fraction of scale, this is for Y axis match between img and screen
heightScaler = 0.4; % 80% height of the screen y axis
imgWidth = 1920 * heightScaler;
imgHeight = imgWidth * orig_width / orig_height;
theRect = [0 0 imgWidth imgHeight];
dstRect_left = CenterRectOnPointd(theRect, screenXpixels/4, screenYpixels/2);
dstRect_right = CenterRectOnPointd(theRect, (screenXpixels/4)*3, screenYpixels/2);

%%

% Load example image. Instruction and question are the same size
theImg=imread(fullfile(pwd,'\C01_TimeDisImg\Img_first\HB_1cut_1_1.png')); % 960 x 540
% Get the size of the image
[orig_width, orig_height, ~]=size(theImg);
% Define the fraction of scale, this is for Y axis match between img and screen
heightScaler = 0.4; % 80% height of the screen y axis
imgWidth = screenXpixels * heightScaler;
imgHeight = imgWidth * orig_width / orig_height;
theRect = [0 0 imgWidth imgHeight];
dstRect_left = CenterRectOnPointd(theRect, screenXpixels/4, screenYpixels/2);
dstRect_right = CenterRectOnPointd(theRect, (screenXpixels/4)*3, screenYpixels/2);