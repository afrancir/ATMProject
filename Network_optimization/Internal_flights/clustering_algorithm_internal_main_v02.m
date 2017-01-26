clc
clear all
close all
tic

% Removing folder NETWORK, if present
if isequal(exist(strcat(pwd,'/NETWORK'),'dir'),7) % 7 = directory
    rmdir(strcat(pwd,'/NETWORK'),'s');
else
end

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

% Creating a list of all sectors
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

fid = fopen(strcat(pwd,'/NETWORK/SECTORS.txt'),'w');
for i=1:numel(sec_IDs)
    fprintf(fid,'%s \n',sec_names{i});
end
fclose(fid);

fid     = fopen(strcat(pwd,'/NETWORK/SECTORS.txt'), 'rt');
SECTORS = textscan(fid,'%s');
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loading all the airports inside the planning domain %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
all_airports = ASPM77.ASPM77(:,1);
all_centers  = ASPM77.ASPM77(:,2);
match_matrix = zeros(length(all_airports),length(CENTERS{1,1}(:)));
for i=1:length(CENTERS{1,1}(:))
    AIRPORTS          = strcmp(CENTERS{1,1}(i),all_centers);
    match_matrix(:,i) = AIRPORTS;
end
% Find all rows characterized by a 1. These are the airports that will be
% considered in the planning domain
[rows,~]=find(match_matrix);
% Checking if the folder containing the airports already exists, and
% creating such folder if teh answer is not
if isequal(exist(strcat(pwd,'/NETWORK/AIRPORTS'),'dir'),7) % 7 = directory
else
    mkdir(strcat(pwd,'/NETWORK/AIRPORTS'));
end
% Create a .txt file with all the airports. Each row corresponds to a
% different airport and has the following entries: airport ID (floating)
% airport name (string), center name (string), latitude (floating), 
% longitude (floating)
fid = fopen(strcat(pwd,'/NETWORK/AIRPORTS/AIRPORTS.txt'),'w');
for i=1:numel(rows)
    fprintf(fid,'%.2f %s %s %.2f %.2f \n',rows(i),ASPM77.ASPM77{rows(i),1},ASPM77.ASPM77{rows(i),2},ASPM77.ASPM77{rows(i),3},ASPM77.ASPM77{rows(i),4});
end
fclose(fid);
% Create 3 additional .txt files. One with the ID-lat-lon of only the
% airports inside the planning domain, one with the ID-lat-lon of all
% airports in the NAS, one with airport ID (floating)
% airport name (string), center name (string), latitude (floating), 
% longitude (floating) of all airports in the NAS
fid = fopen(strcat(pwd,'/NETWORK/AIRPORTS/AIRPORTS_coord.txt'),'w');
for i=1:numel(rows)
    fprintf(fid,'%.2f %.2f %.2f \n',rows(i),ASPM77.ASPM77{rows(i),3},ASPM77.ASPM77{rows(i),4});
end
fclose(fid);
fid = fopen(strcat(pwd,'/NETWORK/AIRPORTS/ALL_AIRPORTS_coord.txt'),'w');
for i=1:numel(all_airports)
    fprintf(fid,'%.2f %.2f %.2f \n',i,ASPM77.ASPM77{i,3},ASPM77.ASPM77{i,4});
end
fclose(fid);
fid = fopen(strcat(pwd,'/NETWORK/AIRPORTS/ALL_AIRPORTS.txt'),'w');
for i=1:numel(all_airports)
    fprintf(fid,'%.2f %s %s %.2f %.2f \n',i,ASPM77.ASPM77{i,1},ASPM77.ASPM77{i,2},ASPM77.ASPM77{i,3},ASPM77.ASPM77{i,4});
end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% For each of the airports considered, store all flight %%%
%%% trajectories that correspond to internal flights      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
internal_flights          = cell(numel(rows),numel(rows));
internal_flights_info     = cell(numel(rows),numel(rows));
cont                      = ones(numel(rows),numel(rows));

time_threshold  = 75;    % [s]
space_threshold = 35;    % [km]: distance covered by an aircraft flying at 550 kts in 2 mins
R_e             = 6371;  % Earth's radius [km]
h               = 10;    % altitude [km]
radius          = R_e+h; % radius of the sphere where we project the points

% Spanning the different files containing flights
for i=1:N_files
    flight_trajectories = load(strcat(pwd,'/TRX_PROCESSED/',listing(i).name));
    flights             = fieldnames(flight_trajectories.ACID_info);
    N_flights           = length(flights);
    % Spanning all the different flights stored in the i-th file
    for j=1:N_flights
        current_flight = flight_trajectories.ACID_info.(flights{j,1});
        if ~isempty(find(strcmp(all_airports,current_flight.origin),1)) && ...
            ismember(find(strcmp(all_airports,current_flight.origin)),rows) && ...
           ~isempty(find(strcmp(all_airports,current_flight.destination),1)) && ...
            ismember(find(strcmp(all_airports,current_flight.destination)),rows)
            
            % ID of the origin airport
            [~,idx_or] = ismember(find(strcmp(all_airports,current_flight.origin)),rows);
            % ID of destination airport
            [~,idx_dest] = ismember(find(strcmp(all_airports,current_flight.destination)),rows);
            % For each time-stamp of the flight, translate the string
            % characterizing the Center into a numerical ID. 1
            % corresponds to the first Center that was specified, and
            % so on so forth. If some of the entries are zero, it means
            % the flight exits the planning domain for a portion of its
            % trajectory
            
            % Checking the Center of reference for each data point of the
            % current flight
            center_mat = zeros(length(current_flight.lat_lon(:,1)),length(CENTERS{1,1}(:)));
            for jj=1:length(CENTERS{1,1}(:))
                center_jj = CENTERS{1,1}(jj);
                center_mat(:,jj) = jj*strcmp(center_jj,current_flight.CurrentCenter);
            end
            % Summing columns to get a vector. Each row will display at
            % most a 1 somewhere, and zeros everywhere else
            center_vec  = sum(center_mat,2);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            sector_mat = zeros(length(current_flight.lat_lon(:,1)),numel(SECTORS{1,1}));
            for jj=1:numel(SECTORS{1,1})
                sector_jj      = SECTORS{1,1}(jj);
                sector_mat(:,jj) = jj*strcmp(sector_jj,current_flight.CurrentSector);
            end
            sector_vec  = sum(sector_mat,2);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            this_flight = horzcat(cont(idx_or,idx_dest)*ones(length(current_flight.lat_lon(:,1)),1),...
                rows(idx_or)*ones(length(current_flight.lat_lon(:,1)),1),...
                rows(idx_dest)*ones(length(current_flight.lat_lon(:,1)),1),...
                current_flight.lat_lon,...
                current_flight.GroundSpeed,current_flight.Heading,current_flight.FlightLevel,...
                current_flight.FiledFL,current_flight.date_vec_UTC,center_vec,sector_vec);
            internal_flights{idx_or,idx_dest} = ...
                vertcat(internal_flights{idx_or,idx_dest},this_flight);
            dummy2      = [{current_flight.ACID_name}, {current_flight.AC_type},...
                {current_flight.origin}, {current_flight.destination},...
                {current_flight.TakeOffDate_UTC}, {current_flight.TakeOffTime_UTC}];
            internal_flights_info{idx_or,idx_dest} = vertcat(internal_flights_info{idx_or,idx_dest},dummy2);
            cont(idx_or,idx_dest)        = cont(idx_or,idx_dest)+1;
        else
        end
    end
end

for i=1:numel(rows)
    if isequal(exist(strcat(pwd,'/NETWORK/INT_traj_info/',num2str(rows(i))),'dir'),7) % 7 = directory
    else
        mkdir(strcat(pwd,'/NETWORK/INT_traj_info/',num2str(rows(i))));
    end
    if isequal(exist(strcat(pwd,'/NETWORK/INT_flight_info/',num2str(rows(i))),'dir'),7) % 7 = directory
    else
        mkdir(strcat(pwd,'/NETWORK/INT_flight_info/',num2str(rows(i))));
    end
    for j=1:numel(rows)
        if ~isempty(internal_flights{i,j})
            dlmwrite(strcat(pwd,'/NETWORK/INT_traj_info/',num2str(rows(i)),'/',num2str(rows(i)),'_',num2str(rows(j)),'_int.txt'),internal_flights{i,j},' ');
            fid = fopen(strcat(pwd,'/NETWORK/INT_flight_info/',num2str(rows(i)),'/',num2str(rows(i)),'_',num2str(rows(j)),'_int.txt'),'w');
            for k=1:length(internal_flights_info{i,j}(:,1))
                fprintf(fid,'%f %s %s %s %s %.0f %.0f %.0f %.0f %.0f %.0f\n',k,internal_flights_info{i,j}{k,1},internal_flights_info{i,j}{k,2},...
                    internal_flights_info{i,j}{k,3},internal_flights_info{i,j}{k,4},...
                    internal_flights_info{i,j}{k,5}(1),internal_flights_info{i,j}{k,5}(2),internal_flights_info{i,j}{k,5}(3),...
                    internal_flights_info{i,j}{k,6}(1),internal_flights_info{i,j}{k,6}(2),internal_flights_info{i,j}{k,6}(3));
            end
            fclose(fid);
        else
        end
    end

end

% Filtering phase: of all the flight that have been stored, keep only
% the ones that satisfy some specific requirements
internal_flights_filtered = cell(numel(rows),numel(rows));
check_flight_filt         = cell(numel(rows),numel(rows));

% Loading airports coordinates
airports                 = load(strcat(pwd,'/NETWORK/AIRPORTS/AIRPORTS_coord.txt'));
% [rad] max. distance with respect to the origin airport. It corresponds to
% roughly 40 km
max_dist_to_airport      = 0.0065; 

% Focusing on the i-th row of cell-array internal_flights
for i=1:numel(rows)
    % Focusing on the j-th column of cell-array internal_flights
    for j=1:numel(rows)
        % at least a flight has been found
        if ~isempty(internal_flights{i,j})
            dummy        = [];
            check_flight = zeros(internal_flights{i,j}(end,1),1);
            for k=1:numel(check_flight)
                % find all entries associated with the j-th flight
                idx    = find(internal_flights{i,j}(:,1)==k);
                % if there are less than 10 entries overall,
                % do not consider the flight
                if numel(idx)<10
                else
                    flight = internal_flights{i,j}(idx,:);
                    % now, store only the entries inside the planning domain
                    idx2          = find(flight(:,end-1)~=0);
                    % if there are less than 10 entries in the planning
                    % domain, do not consider the flight
                    if numel(idx2)<10
                    else
                        flight_domain = flight(idx2,:);
                        % for these entries, compute the time-gap between consecutive
                        % entries
                        time_difference = etime(flight_domain(2:end,end-7:end-2),flight_domain(1:end-1,end-7:end-2));
                        % if the maximum time-distance between consecutive entries
                        % if less than a specified threshold, the flight is kept
                        if max(time_difference)<time_threshold
                            % Compute distance between consecutive entries to make
                            % sure there are no spatial "gaps" that make no
                            % physical sense
                            dist_consecutive_entries = zeros(numel(idx2)-1,1);
                            for kk=1:numel(dist_consecutive_entries)
                                dist_consecutive_entries(kk,1) = Haversine(flight_domain(kk+1,4:5),flight_domain(kk,4:5),1)*radius;
                            end
                            if any(dist_consecutive_entries>space_threshold)
                            else
                                % Check if the trajectory starts sufficiently
                                % close to the origin airport and ends 
                                % sufficiently close to the destination 
                                % airport
                                idx_or   = find(airports(:,1)==flight(1,2));
                                idx_dest = find(airports(:,1)==flight(1,3));
                                if (Haversine(airports(idx_or,2:3),flight(idx2(1),4:5),1)<=max_dist_to_airport && ...
                                    Haversine(airports(idx_dest,2:3),flight(idx2(end),4:5),1)<=max_dist_to_airport)    
                                    check_flight(k) = 1;
                                    dummy           = [dummy;flight_domain];
                                else
                                end
                            end
                        else
                        end
                    end
                end
            end
            internal_flights_filtered{i,j} = dummy;
            check_flight_filt{i,j}         = check_flight;
        else
        end
    end
end

for i=1:numel(rows)
    if isequal(exist(strcat(pwd,'/NETWORK/INT_traj_filt_info/',num2str(rows(i))),'dir'),7) % 7 = directory
    else
        mkdir(strcat(pwd,'/NETWORK/INT_traj_filt_info/',num2str(rows(i))));
    end
    if isequal(exist(strcat(pwd,'/NETWORK/INT_flight_filt_info/',num2str(rows(i))),'dir'),7) % 7 = directory
    else
        mkdir(strcat(pwd,'/NETWORK/INT_flight_filt_info/',num2str(rows(i))));
    end
    for j=1:numel(rows)
        if ~isempty(internal_flights_filtered{i,j})
            dlmwrite(strcat(pwd,'/NETWORK/INT_traj_filt_info/',num2str(rows(i)),'/',num2str(rows(i)),'_',num2str(rows(j)),'_int.txt'),internal_flights_filtered{i,j},' ');
            fid = fopen(strcat(pwd,'/NETWORK/INT_flight_filt_info/',num2str(rows(i)),'/',num2str(rows(i)),'_',num2str(rows(j)),'_int.txt'),'w');
            for k=1:length(check_flight_filt{i,j}(:,1))
                if check_flight_filt{i,j}(k,1)
                    fprintf(fid,'%f %s %s %s %s %.0f %.0f %.0f %.0f %.0f %.0f\n',k,internal_flights_info{i,j}{k,1},internal_flights_info{i,j}{k,2},...
                    internal_flights_info{i,j}{k,3},internal_flights_info{i,j}{k,4},...
                    internal_flights_info{i,j}{k,5}(1),internal_flights_info{i,j}{k,5}(2),internal_flights_info{i,j}{k,5}(3),...
                    internal_flights_info{i,j}{k,6}(1),internal_flights_info{i,j}{k,6}(2),internal_flights_info{i,j}{k,6}(3));
                else
                end
            end
            fclose(fid);
        else
        end
    end
    
end

time = toc;