function[TTLdescription]=TTLinterpret(TTLnumber)

%% function[TTLdescription]=TTLinterpret(TTLnumber)

%% Create table for TTL value key
%TTLvalue={'61','11','1','2','3','4','60','0'};
TTLvalue={'61','11','1','2','3','111','112','113','114','115','116','60','0'};

%TTL_ID={'TaskOnset','FixCross','ClipOnset','ClipOffset','Probe','Response','TaskOffset','AfterTTLdelay'};
TTL_ID={'TaskOnset','FixCross','ClipOnset','ClipOffset','Probe','Response','Response',...
    'Response','Response','Response','Response','TaskOffset','AfterTTLdelay'};

%% Create key table
keyTable=table(transpose(TTL_ID),transpose(TTLvalue),'VariableNames',{'TTLI','TTLN'});

%% Create empty array to preallocate
TTLdescription=cell(size(TTLnumber));

%% Loop through TTLs
for i=1:length(TTLnumber)
    tempN=TTLnumber{i};

    %% Find matching ID from keyTable if it exists
    if ~ismember(tempN,keyTable.TTLN)
        TTLdescription{i} = 'NaN';
    else
        TTLlabel=keyTable.TTLI(ismember(keyTable.TTLN,tempN));
        TTLdescription{i}=TTLlabel{1};
    end
end

end