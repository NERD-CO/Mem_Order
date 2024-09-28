function [cleanPOS] = cleanUPpos_MO(posIN , type , xbounds , ybounds)

% long values
% posL1 = posT2(posT2(:,1) > 1,:);
% posL2 = posL1(posL1(:,2) > 1,:);

switch type
    case 1
        cleanPOS = double(posIN);
        cleanPOS(cleanPOS > 6000) = nan;
    case 2
        inframeTEST = zeros(height(posIN),2,'logical');
        tmpDOUBLE = double(posIN);
        inframeTEST(:,1) = tmpDOUBLE(:,1) > 0 & tmpDOUBLE(:,1) < 1920;
        inframeTEST(:,2) = tmpDOUBLE(:,2) > 0 & tmpDOUBLE(:,2) < 1320;
        cleanPOS = all(inframeTEST,2);

        testPOS = posIN(cleanPOS,:);
        plot(testPOS(:,1),testPOS(:,2));
    case 3
        close all
        inframeTEST = zeros(height(posIN),2,'logical');
        tmpDOUBLE = double(posIN);
        inframeTEST(:,1) = tmpDOUBLE(:,1) > xbounds(1) &...
            tmpDOUBLE(:,1) < xbounds(2);
        inframeTEST(:,2) = tmpDOUBLE(:,2) > (ybounds(1)) &...
            tmpDOUBLE(:,2) < (ybounds(2));
        cleanPOS = all(inframeTEST,2);

        tiledlayout(1,2)
        testPOS = posIN(cleanPOS,:);
        % plot(posIN(:,2),posIN(:,1));ylim([0 1080]);xlim([0 1920])
        nexttile
        plot(testPOS(:,1),testPOS(:,2))
        xlim([0 1920])
        ylim([0 1320])
        xline(xbounds)
        yline(ybounds)

        rawPOS = double(posIN);
        rawPOS(rawPOS > 6000) = nan;
        nexttile
        plot(rawPOS(:,1),rawPOS(:,2),'r')
        xlim([0 1920])
        ylim([0 1320])
        xline(xbounds)
        yline(ybounds)
        pause(0.5)



end

end