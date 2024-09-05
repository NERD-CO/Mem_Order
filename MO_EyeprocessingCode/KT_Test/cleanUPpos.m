function [cleanPOS] = cleanUPpos(posIN)

%% function [cleanPOS] = cleanUPpos(posIN)

%%
posS1 = smoothdata(posIN(:,1),'gaussian',40);
posS2 = smoothdata(posIN(:,2),'gaussian',40);
cleanPOS = [posS1 , posS2];

end