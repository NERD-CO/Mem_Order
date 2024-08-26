%% Clear
clear;
clc;

%% Scratch
file3 = nwbRead('Z:\Tyner_K_Projects\MEM_ORDER\DataFolder\MW19\NWBProcessing\NWB_Data\MW19_Session_3_filter.nwb');
file4 = nwbRead('Z:\Tyner_K_Projects\MEM_ORDER\DataFolder\MW19\NWBProcessing\NWB_Data\MW19_Session_4_filter.nwb');
file12 = nwbRead('Z:\Tyner_K_Projects\MEM_ORDER\DataFolder\MW19\NWBProcessing\NWB_Data\MW19_Session_12_filter.nwb');
file14 = nwbRead('Z:\Tyner_K_Projects\MEM_ORDER\DataFolder\MW19\NWBProcessing\NWB_Data\MW19_Session_14_filter.nwb');
file15 = nwbRead('Z:\Tyner_K_Projects\MEM_ORDER\DataFolder\MW19\NWBProcessing\NWB_Data\MW19_Session_15_filter.nwb');

%% Grab TTLs
file3_et = file3.acquisition.get('events').timestamps.load();
file3_ttl = cellstr(file3.acquisition.get('events').data.load());
I = contains(file3_ttl,'TTL');
file3_et = file3_et(I);
file3_ttl = file3_ttl(I);
file3_ttl = extractBetween(file3_ttl,'(',')');
file3_ttl = cellfun(@(x) hex2dec(x),file3_ttl,'UniformOutput',true);

file4_et = file4.acquisition.get('events').timestamps.load();
file4_ttl = cellstr(file4.acquisition.get('events').data.load());
I = contains(file4_ttl,'TTL');
file4_et = file4_et(I);
file4_ttl = file4_ttl(I);
file4_ttl = extractBetween(file4_ttl,'(',')');
file4_ttl = cellfun(@(x) hex2dec(x),file4_ttl,'UniformOutput',true);

file12_et = file12.acquisition.get('events').timestamps.load();
file12_ttl = cellstr(file12.acquisition.get('events').data.load());
I = contains(file12_ttl,'TTL');
file12_et = file12_et(I);
file12_ttl = file12_ttl(I);
file12_ttl = extractBetween(file12_ttl,'(',')');
file12_ttl = cellfun(@(x) hex2dec(x),file12_ttl,'UniformOutput',true);

file14_et = file14.acquisition.get('events').timestamps.load();
file14_ttl = cellstr(file14.acquisition.get('events').data.load());
I = contains(file14_ttl,'TTL');
file14_et = file14_et(I);
file14_ttl = file14_ttl(I);
file14_ttl = extractBetween(file14_ttl,'(',')');
file14_ttl = cellfun(@(x) hex2dec(x),file14_ttl,'UniformOutput',true);

file15_et = file15.acquisition.get('events').timestamps.load();
file15_ttl = cellstr(file15.acquisition.get('events').data.load());
I = contains(file15_ttl,'TTL');
file15_et = file15_et(I);
file15_ttl = file15_ttl(I);
file15_ttl = extractBetween(file15_ttl,'(',')');
file15_ttl = cellfun(@(x) hex2dec(x),file15_ttl,'UniformOutput',true);

diff1 = (file4_et(1,1) - file3_et(end,1))/(1e6);
diff2 = (file12_et(1,1) - file4_et(end,1))/(1e6);
diff3 = (file14_et(1,1) - file12_et(end,1))/(1e6);
diff4 = (file15_et(1,1) - file14_et(end,1))/(1e6);