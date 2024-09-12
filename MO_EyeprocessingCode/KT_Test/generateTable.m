function outTable = generateTable(respMat,task)

%% function outTable = generateTable(respmat)

%% Grab cross times
CrossStart = [respMat.CrossStart]';
CrossEnd = [respMat.CrossEnd]';

%% Grabe Clip/Frame times
if isfield(respMat,'ClipStart')
    FrameStart = [respMat.ClipStart]';
elseif isfield(respMat,'FrameStart')
    FrameStart = [respMat.FrameStart]';
end

if isfield(respMat,'ClipEnd')
    FrameEnd = [respMat.ClipEnd]';
elseif isfield(respMat,'FrameEnd')
    FrameEnd = [respMat.FrameEnd]';
end

%% Grab Question times
QuesStart = [respMat.QuesStart]';
QuesEnd = [respMat.respTime]';

%% Grab respValue
respValue = NaN(length(respMat),1);
for ii = 1:length(respMat)
    if ~isempty(respMat(ii).respValue)
        respValue(ii,1) = respMat(ii).respValue;
    else
        respValue(ii,1) = NaN;
    end
end

%% Frame Name
if isfield(respMat,'ClipName')
    FrameName = {respMat.ClipName}';
elseif isfield(respMat,'FrameName')
    FrameName = {respMat.FrameName}';
end

%% Question Name
if isfield(respMat,'QuesName')
    QuesName = {respMat.QuesName}';
else
    QuesName = repmat({NaN},length(CrossStart),1);
end

%% Frame Order
if isfield(respMat,'FrameOrder')
    FrameOrder = {respMat.FrameOrder}';
else
    FrameOrder = repmat({NaN},length(CrossStart),1);
end

%% Pt responses
switch task
    case 'encoding'
        PtCorrect = NaN(length(CrossStart),1);
        RealAnswer = repmat({NaN},length(CrossStart),1);
    case 'sceneRecog'
        %% Identify correct responses
        % 1 = correct response, 0 = incorrect response
        I = contains(FrameName,'foil');
        if any(I)
            RealAnswer = ~contains(FrameName,'foil');
        else
            RealAnswer = repmat({NaN},length(CrossStart),1);
        end

        %% Identify Pt correct
        PtCorrect = NaN(length(RealAnswer),1);
        if ~iscell(RealAnswer)
            for ii = 1:length(RealAnswer)
                if RealAnswer(ii,1) == 1 && ((respValue(ii,1)== 111) || (respValue(ii,1) == 112) || (respValue(ii,1) == 113))
                    PtCorrect(ii,1) = 1;
                elseif RealAnswer(ii,1) == 1 && ((respValue(ii,1)== 114) || (respValue(ii,1) == 115) || (respValue(ii,1) == 116))
                    PtCorrect(ii,1) = 0;
                elseif RealAnswer(ii,1) == 0 && ((respValue(ii,1)== 111) || (respValue(ii,1) == 112) || (respValue(ii,1) == 113))
                    PtCorrect(ii,1) = 0;
                elseif RealAnswer(ii,1) == 0 && ((respValue(ii,1)== 114) || (respValue(ii,1) == 115) || (respValue(ii,1) == 116))
                    PtCorrect(ii,1) = 1;
                end
            end
        end
    case 'timeDiscrim'
        %% Identify correct responses
        I = cellfun(@(c) any(isnan(c(:))), FrameOrder);
        RealAnswer = cell(length(FrameOrder),1);
        if ~any(I)
            for ii = 1:length(FrameOrder)
                order = FrameOrder{ii,1};
                if order(1,1) == 2
                    RealAnswer{ii,1} = 'R';
                else
                    RealAnswer{ii,1} = 'L';
                end
            end
        end

        %% Identify Pt correct
        PtCorrect = NaN(length(FrameOrder),1);
        for ii = 1:length(FrameOrder)
            order = FrameOrder{1,1};
            if order(1,1) == 2 && ((respValue(ii,1) == 111) || (respValue(ii,1) == 112) || (respValue(ii,1) == 113))
                PtCorrect(ii,1) = 1;
            elseif order(1,1) == 1 && ((respValue(ii,1) == 111) || (respValue(ii,1) == 112) || (respValue(ii,1) == 113))
                PtCorrect(ii,1) = 0;
            elseif order(1,1) == 2 && ((respValue(ii,1) == 114) || (respValue(ii,1) == 115) || (respValue(ii,1) == 116))
                PtCorrect(ii,1) = 0;
            elseif order(1,1) == 1 && ((respValue(ii,1) == 114) || (respValue(ii,1) == 115) || (respValue(ii,1) == 116))
                PtCorrect(ii,1) = 1;
            end
        end
end

%% Get the boundary type
Boundary = extractBefore(FrameName,'_');

%% Build Table
outTable = table(CrossStart,CrossEnd,FrameStart,FrameEnd,FrameName,FrameOrder,QuesStart,QuesEnd,QuesName,respValue,RealAnswer,PtCorrect,Boundary);

end