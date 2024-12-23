function [] = summarizeCSVs_MO(CSV_FOLDER)

% CD to CSV folder
% Load CSV files in aggregate
% Tabluate the number of unique hemispheres/regions/contact pairs per
% subject
close all
cd(CSV_FOLDER)

csvDIR = dir('*.xlsx');
csvDIR2 = {csvDIR.name};

allTAB = table;
for cii = 1:length(csvDIR2)
    tmpTab = readtable(csvDIR2{cii});
    subjectID = csvDIR2{cii}(1:4);
    tmpTab.subject = repmat({subjectID},height(tmpTab),1);
    tmpTab = movevars(tmpTab,"subject","Before","ch_name");
    allTAB = [allTAB ; tmpTab];
end

% Check all
itpcPresent = any(logical(table2array(allTAB(:,...
    ["first","second","third","fourth"]))),2);

subRegITPC = allTAB(itpcPresent,:);

% Number of positive ITPC subjects
numITPCsubs = numel(unique(subRegITPC.subject));
fracITPCsubs = numITPCsubs/8;

% Number / Fraction of regions / contacts per subject
uniIDSsub = unique(subRegITPC.subject);
regionSummary = zeros(numITPCsubs,2);
contactSummary = zeros(numITPCsubs,2);
for nii = 1:numITPCsubs

    tmpSubid = uniIDSsub(nii);

    itpcTab = subRegITPC(matches(subRegITPC.subject,tmpSubid),:);
    totTab = allTAB(matches(allTAB.subject,tmpSubid),:);

    contactSummary(nii,1) = height(itpcTab);
    contactSummary(nii,2) = height(itpcTab)/height(totTab);

    itpcREGs = numel(unique(itpcTab.ch_longName));
    totREGs = numel(unique(totTab.ch_longName));

    regionSummary(nii,1) = itpcREGs;
    regionSummary(nii,2) = itpcREGs/totREGs;

end

newcolors = [0.83 0.14 0.14
             1.00 0.54 0.00
             0.47 0.25 0.80
             0.25 0.80 0.54
             0 1 0 
             1 0 0
             0.5 0 1
             0 0 0
             0 0 1];
         
% Donut summaries
figure;
tiledlayout(1,2)
nexttile
d1 = donutchart(regionSummary(:,1));
d1.Names = uniIDSsub;
% d1.Labels = [transpose(d1.Names) ,'[', string(regionSummary(:,1)),']'];
title('Brain Area expression of ITPC')
colororder(newcolors)

nexttile
d2 = donutchart(contactSummary(:,1));
d2.Names = uniIDSsub;
% d1.Labels = [transpose(d1.Names) ,'[', string(regionSummary(:,1)),']'];
title('Contact # expression of ITPC')
colororder(newcolors)

set(gcf,'Position',[218 254 1171 1032])

% Plot overall region IDs
uniREGS = unique(subRegITPC.ch_longName);

% number of contacts per reg (across subjects)
% number of subjects sharing reg

regAcross = zeros(numel(uniREGS),4);
for uuir = 1:numel(uniREGS)

    tmpREGi = uniREGS{uuir};

    regAll = subRegITPC(matches(subRegITPC.ch_longName,tmpREGi),:);

    regAcross(uuir,1) = height(regAll);
    regAcross(uuir,2) = height(regAll)/height(subRegITPC);

    tmpSUBci = numel(unique(regAll.subject));

    regAcross(uuir,3) = tmpSUBci;
    regAcross(uuir,4) = tmpSUBci/numITPCsubs;
end

% Get group by summary of region by subject
% Limit by greater than 2 subjects sharing region
regAcross2 = regAcross(regAcross(:,3) > 1,:);
regAcrossIDs2 = uniREGS(regAcross(:,3) > 1);


figure;
tiledlayout(1,2)
nexttile
d3 = donutchart(regAcross2(:,1));
d3.Names = regAcrossIDs2;
% d1.Labels = [transpose(d1.Names) ,'[', string(regionSummary(:,1)),']'];
title('Fraction of contacts per region')
colororder(newcolors)

nexttile
d4 = donutchart(regAcross2(:,3));
d4.Names = regAcrossIDs2;
d4.LabelStyle = "namedata";
% d1.Labels = [transpose(d1.Names) ,'[', string(regionSummary(:,1)),']'];
title('Subject # sharing region')
colororder(newcolors)

set(gcf,'Position',[218 214 1458 1072])






end