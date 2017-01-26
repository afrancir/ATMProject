clc
clear all
close all
tic

ASPM77              = load(strcat(pwd,'/AUX_DATA/ASPM77.mat'));
listing             = dir(strcat(pwd,'/TRX_PROCESSED/*.mat'));
N_files             = 14;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loading the Centers that form the planning domain %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid     = fopen(strcat(pwd,'/AUX_DATA/PLANNING_DOMAIN.txt'), 'rt');
CENTERS = textscan(fid,'%s');
fclose(fid);

% Loading names of all Centers in the NAS
all_info     = load(strcat(pwd,'/AUX_DATA/waypoint_info_v01.mat'));
Center_names = fieldnames(all_info.center_bound);

if isequal(exist(strcat(pwd,'/NETWORK/Center_boundaries'),'dir'),7) % 7 = directory
else
    mkdir(strcat(pwd,'/NETWORK/Center_boundaries'));
end

for i=1:length(CENTERS{1,1}(:,1))
    current_Center = CENTERS{1,1}(i,1);
    idx            = find(strcmp(Center_names,current_Center));
    current_Center_latlon = all_info.center_bound.(current_Center{1,1});
    dlmwrite(strcat(pwd,'/NETWORK/Center_boundaries/',current_Center{1,1},'.txt'),current_Center_latlon,' ');
end

names     = fieldnames(all_info.sector_info);
cont      = 1;
sec_names = [];
sec_IDs   = [];

for i=1:length(CENTERS{1,1}(:,1))
    current_Center = CENTERS{1,1}(i,1);
    index          = strmatch(current_Center,names);
    this_sec_ID    = cont:cont+numel(index)-1;
    cont           = this_sec_ID(end)+1;
    sec_names      = [sec_names;names(index)];
    sec_IDs        = [sec_IDs;this_sec_ID'];
end

fid = fopen(strcat(pwd,'/SECTORS.txt'),'w');
for i=1:numel(sec_IDs)
    fprintf(fid,'%s \n',sec_names{i});
end
fclose(fid);

