% Written by Alessandro Bombelli, 22th January 2017
% Code that analyzes recorded trajectories within the planning domain
% selected, and stores hourly departures/arrivals for the airports of
% interest to compute Gilbo envelopes

clc
clear all
close all

ASPM77              = load(strcat(pwd,'/Internal_flights/AUX_DATA/ASPM77.mat'));
listing             = dir(strcat(pwd,'/Internal_flights/TRX_PROCESSED/*.mat'));
N_files             = 1;


all_airports = ASPM77.ASPM77(:,1);

% Airports within planning domain. Each airport is characterized by an
% array with 3 components, i.e., airport_ID, latitude, longitude
airports_pd    = load(strcat(pwd,'/Internal_flights/NETWORK/AIRPORTS/AIRPORTS_coord.txt'));
airports_pd_ID = airports_pd(:,1);

departures     = cell(numel(airports_pd_ID),1);
arrivals       = cell(numel(airports_pd_ID),1);
cont_dep       = ones(numel(airports_pd_ID),1);
cont_arr       = ones(numel(airports_pd_ID),1);


for i=1:N_files
    flight_trajectories = load(strcat(pwd,'/Internal_flights/TRX_PROCESSED/',listing(i).name));
    flights             = fieldnames(flight_trajectories.ACID_info);
    N_flights           = length(flights);
    for j=1:N_flights
        current_flight = flight_trajectories.ACID_info.(flights{j,1});
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% ID of the origin airport is consistent with one of the %%%
        %%% airports we are considering                            %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isempty(find(strcmp(all_airports,current_flight.origin),1)) && ...
                ismember(find(strcmp(all_airports,current_flight.origin)),airports_pd_ID)
            
            % Detect which airport we are considering
            [~,idx_or] = ismember(find(strcmp(all_airports,current_flight.origin)),airports_pd_ID);
            % Lat-Lon of all entries associated with this flight
            Lat_Lon     = current_flight.lat_lon;
            % Lat-Lon of first entry
            Lat_Lon_dep = Lat_Lon(1,:);
            % Lat-Lon of origin airport
            Lat_Lon_origin_airport = airports_pd(airports_pd(:,1)==airports_pd_ID(idx_or),2:3);
            % Check if first recorded entry is sufficiently close to the
            % origin airport
            if Haversine(Lat_Lon_dep,Lat_Lon_origin_airport,1)<=40/6378
                departures{idx_or,1} = vertcat(departures{idx_or,1},horzcat(cont_dep(idx_or),...
                    current_flight.TakeOffDate_UTC,current_flight.TakeOffTime_UTC));
                cont_dep(idx_or)         = cont_dep(idx_or)+1;
            else
            end
        else
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% ID of the destination airport is consistent with one of the %%%
        %%% airports we are considering                                 %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isempty(find(strcmp(all_airports,current_flight.destination),1)) && ...
                ismember(find(strcmp(all_airports,current_flight.destination)),airports_pd_ID)
            
            % Detect which airport we are considering
            [~,idx_dest] = ismember(find(strcmp(all_airports,current_flight.destination)),airports_pd_ID);
            % Lat-Lon of all entries associated with this flight
            Lat_Lon     = current_flight.lat_lon;
            % Lat-Lon of last entry
            Lat_Lon_dep = Lat_Lon(end,:);
            % Lat-Lon of destination airport
            Lat_Lon_destination_airport = airports_pd(airports_pd(:,1)==airports_pd_ID(idx_dest),2:3);
            % Check if last recorded entry is sufficiently close to the
            % destination airport
            if Haversine(Lat_Lon_dep,Lat_Lon_destination_airport,1)<=40/6378
                arrivals{idx_dest,1} = vertcat(arrivals{idx_dest,1},horzcat(cont_arr(idx_dest),...
                    current_flight.date_vec_UTC(end,:)));
                cont_arr(idx_dest)         = cont_arr(idx_dest)+1;
            else
            end
        else
        end
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Creating folder DEPARTURES_Gilbo %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isequal(exist(strcat(pwd,'/DEPARTURES_Gilbo'),'dir'),7) % 7 = directory
    mkdir(strcat(pwd,'/DEPARTURES_Gilbo'));
else
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Creating folder ARRIVALS_Gilbo %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isequal(exist(strcat(pwd,'/ARRIVALS_Gilbo'),'dir'),7) % 7 = directory
    mkdir(strcat(pwd,'/ARRIVALS_Gilbo'));
else
end

for i=1:numel(airports_pd_ID)
    if ~isempty(departures{i,1})
        dlmwrite(strcat(pwd,'/DEPARTURES_Gilbo/',num2str(airports_pd_ID(i)),'_.txt'),departures{i,1},' ');
    else
    end
    if ~isempty(arrivals{i,1})
        dlmwrite(strcat(pwd,'/ARRIVALS_Gilbo/',num2str(airports_pd_ID(i)),'_.txt'),arrivals{i,1},' ');
    else
    end
end

for i=3
    if ~isempty(load(strcat(pwd,'/DEPARTURES_Gilbo/',num2str(airports_pd_ID(i)),'_.txt')))
        Departures_from_this_airport = load(strcat(pwd,'/DEPARTURES_Gilbo/',num2str(airports_pd_ID(i)),'_.txt'));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Determine which days have been processed    %%%
        %%% as it concerns departures from this airport %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        YMD_dep          = unique(Departures_from_this_airport(:,2:4),'rows');
        Dep_mat          = zeros(length(YMD_dep(:,1)),27);
        Dep_mat(:,25:27) = YMD_dep;
        for j=1:numel(Departures_from_this_airport(:,1))
            [~,row_day]   = ismember(Departures_from_this_airport(j,2:4),YMD_dep,'rows');
            if row_day ~= 0
                this_hour                    = Departures_from_this_airport(j,5);
                Dep_mat(row_day,this_hour+1) = Dep_mat(row_day,this_hour+1)+1;
            else
            end
        end
        
        dlmwrite(strcat(pwd,'/DEPARTURES_Gilbo/',num2str(airports_pd_ID(i)),'_DEP.txt'),Dep_mat,' ');
        
    else
    end
    
    if ~isempty(load(strcat(pwd,'/ARRIVALS_Gilbo/',num2str(airports_pd_ID(i)),'_.txt')))
        Arrivals_from_this_airport = load(strcat(pwd,'/ARRIVALS_Gilbo/',num2str(airports_pd_ID(i)),'_.txt'));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Determine which days have been processed  %%%
        %%% as it concerns arrivals from this airport %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        YMD_arr          = unique(Arrivals_from_this_airport(:,2:4),'rows');
        Arr_mat          = zeros(length(YMD_arr(:,1)),27);
        Arr_mat(:,25:27) = YMD_arr;
        for j=1:numel(Arrivals_from_this_airport(:,1))
            [~,row_day]   = ismember(Arrivals_from_this_airport(j,2:4),YMD_arr,'rows');
            if row_day ~= 0
                this_hour                    = Arrivals_from_this_airport(j,5);
                Arr_mat(row_day,this_hour+1) = Arr_mat(row_day,this_hour+1)+1;
            else
            end
        end
        
        dlmwrite(strcat(pwd,'/ARRIVALS_Gilbo/',num2str(airports_pd_ID(i)),'_ARR.txt'),Arr_mat,' ');
        
    else
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Creating folder GILBO_envelopes %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isequal(exist(strcat(pwd,'/GILBO_envelopes'),'dir'),7) % 7 = directory
    mkdir(strcat(pwd,'/GILBO_envelopes'));
else
end

for i=3
    % Check if the airport has both a departure matrix and an arrival
    % matrix
    if ~isempty(load(strcat(pwd,'/DEPARTURES_Gilbo/',num2str(airports_pd_ID(i)),'_DEP.txt'))) && ...
            ~isempty(load(strcat(pwd,'/ARRIVALS_Gilbo/',num2str(airports_pd_ID(i)),'_ARR.txt')))
        
        GILBO           = [];
        
        dep_MAT         = load(strcat(pwd,'/DEPARTURES_Gilbo/',num2str(airports_pd_ID(i)),'_DEP.txt'));
        arr_MAT         = load(strcat(pwd,'/ARRIVALS_Gilbo/',num2str(airports_pd_ID(i)),'_ARR.txt'));
        days_dep        = dep_MAT(:,25:27);
        days_arr        = arr_MAT(:,25:27);
        [~,i_dep,i_arr] = intersect(days_dep,days_arr,'rows');
        
        % We found at least one common day 
        if ~isempty(i_dep)
            for j=1:numel(i_dep)
                this_day_dep = dep_MAT(i_dep(j),1:24);
                this_day_arr = arr_MAT(i_arr(j),1:24);
                this_day     = horzcat(this_day_dep',this_day_arr',(0:23)',repmat(dep_MAT(i_dep(j),25:27),24,1));
                GILBO        = vertcat(GILBO,this_day);
                
                dlmwrite(strcat(pwd,'/GILBO_envelopes/',num2str(airports_pd_ID(i)),'_Gilbo.txt'),GILBO,' ');
                
            end
        else
        end
        
    else
    end
end


for i=3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Loading GILBO Matrix and computing Convex Hull %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ~isempty(load(strcat(pwd,'/GILBO_envelopes/',num2str(airports_pd_ID(i)),'_Gilbo.txt')))
        GILBO_this_airp = load(strcat(pwd,'/GILBO_envelopes/',num2str(airports_pd_ID(i)),'_Gilbo.txt'));
        Dep_Arr         = GILBO_this_airp(:,1:2);
        max_Dep         = max(Dep_Arr(:,1));
        max_Arr         = max(Dep_Arr(:,2));
        aux_point1      = [0 0];
        aux_point2      = [max_Dep 0];
        aux_point3      = [0 max_Arr];
        aux_points      = vertcat(aux_point1,aux_point2,aux_point3);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Plot data and auxiliary points %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure()
        hold on
        plot(Dep_Arr(:,1),Dep_Arr(:,2),'Color','b','Linestyle','none','Linewidth',2,'Marker','s','Markersize',8)
        plot(aux_points(:,1),aux_points(:,2),'Color','k','Linestyle','none','Linewidth',2,'Marker','x','Markersize',10)
        xlabel('Departures [1/h]','Fontname','Avantgarde','Fontsize',14)
        ylabel('Arrivals [1/h]','Fontname','Avantgarde','Fontsize',14)
        grid on
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Determine convex hull %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        all_data    = vertcat(Dep_Arr,aux_points);
        idx_cvxhull = convhull(all_data);
        cvxhull     = all_data(idx_cvxhull,:);
        plot(cvxhull(:,1),cvxhull(:,2),'Color','r','Linestyle','-','Linewidth',2,'Marker','s','Markersize',8)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Determine horizontal boundary on max arrivals %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        idx_zero_dep   = find(cvxhull(:,1)==0);
        max_arrivals   = max(cvxhull(idx_zero_dep,2));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Determine vertical boundary on max departures %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        idx_zero_arr   = find(cvxhull(:,2)==0);
        max_departures = max(cvxhull(idx_zero_arr,1));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Determine sloped constraint %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        idx_both_dep_arr = setdiff((1:length(cvxhull(:,1)))',vertcat(idx_zero_dep,idx_zero_arr));
        Y                = cvxhull(idx_both_dep_arr,2);
        X                = horzcat(ones(numel(idx_both_dep_arr),1),cvxhull(idx_both_dep_arr,1));
        Beta             = X\Y;
        
    else
    end
end



