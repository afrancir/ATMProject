% Written by Alessandro Bombelli, August 19, 2016
% Given a time-step

clc
clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loading the Centers that form the planning domain %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid     = fopen(strcat(pwd,'/AUX_DATA/PLANNING_DOMAIN.txt'), 'rt');
CENTERS = textscan(fid,'%s');
fclose(fid);

AIRPORTS   = load(strcat(pwd,'/NETWORK/AIRPORTS/AIRPORTS_coord.txt'));
airport_ID = AIRPORTS(:,1);

% Checking if the folder exists already or not
if isequal(exist(strcat(pwd,'/NETWORK/INT_AggRoutes_timestep'),'dir'),7) % 7 = directory
else
    mkdir(strcat(pwd,'/NETWORK/INT_AggRoutes_timestep'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Selection of the time-step %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 5;

% Origin airport
for i=1:numel(airport_ID)
    
    % Checking if the folder exists already or not
    if isequal(exist(strcat(pwd,'/NETWORK/INT_AggRoutes_timestep/',num2str(airport_ID(i))),'dir'),7) % 7 = directory
    else
        mkdir(strcat(pwd,'/NETWORK/INT_AggRoutes_timestep/',num2str(airport_ID(i))));
    end
    
    AGG_R_dt     = [];
    cont         = 1;
    
    % Destination airport
    for k=1:numel(airport_ID)
        check_this_OD = isequal(exist(strcat(pwd,'/NETWORK/INT_AggRoutes/',num2str(airport_ID(i)),'/',num2str(airport_ID(i)),'_',num2str(airport_ID(k)),'.txt'),'file'),2);
        
        if check_this_OD
            
            AGG_R        = load(strcat(pwd,'/NETWORK/INT_AggRoutes/',num2str(airport_ID(i)),'/',num2str(airport_ID(i)),'_',num2str(airport_ID(k)),'.txt'));
            N_agg_routes = AGG_R(end,1);
            
            
            for j=1:N_agg_routes
                % find all rows associated with current aggregate route
                idx        = find(AGG_R(:,1)==j);
                this_route = AGG_R(idx,:);
                % get length of the route
                N_r   = numel(idx);
                n_bar = floor((N_r-1)/(2*dt));
                eta1  = 1+n_bar*dt;    % last node starting from the origin node
                eta2  = N_r-n_bar*dt;  % first node starting from the destination node
                eta   = eta2-eta1;     % difference between the two nodes (ranges from 0 to 2*dt-1)
                
                % do nothing in this case
                if eta == 0
                    
                    idx_left      = (1:dt:1+n_bar*dt)';
                    idx_right     = (N_r-n_bar*dt:dt:N_r)';
                    idx_dt        = unique(vertcat(idx_left,idx_right));
                    this_route_dt = this_route(idx_dt,2:end);
                    
                    % eliminate eta1 and eta2 and define an auxiliary point in
                    % between
                    
                elseif eta >= 1 && eta <= floor(dt/2)
                    
                    % if n_bar=0, we have a very short route (probably we are
                    % very close to the boundary). Thus, we do not create an
                    % intermediate node and delete eta1 and eta2, but simply
                    % keep the two extremes
                    if n_bar == 0
                        this_route_dt = this_route([1;N_r],:);
                    else
                        idx_left            = (1:dt:1+(n_bar-1)*dt)';
                        idx_right           = (N_r-(n_bar+1)*dt:dt:N_r)';
                        this_route_dt_left  = this_route(idx_left,:);
                        this_route_dt_right = this_route(idx_right,:);
                        info_eta1           = this_route(eta1,:);
                        info_eta2           = this_route(eta2,:);
                        dummy               = zeros(1,length(this_route(1,:)));
                        dummy(1)            = info_eta1(1);
                        dummy(2)            = info_eta1(2);
                        dummy(3)            = info_eta1(3);
                        [lat,lon]           = gcwaypts(info_eta1(4),info_eta1(5),info_eta2(4),info_eta2(5),3);
                        dummy(4)            = lat(2);
                        dummy(5)            = lon(2);
                        dummy(6)            = info_eta1(6);
                        dummy(7)            = info_eta1(7);
                        dummy(8)            = info_eta1(8);
                        dummy(9)            = info_eta1(9);
                        dummy(10)           = info_eta1(10);
                        dummy(11)           = info_eta1(11);
                        
                        this_route_dt       = vertcat(this_route_dt_left(:,2:end),dummy(2:end),this_route_dt_right(:,2:end));
                    end
                    
                    % do nothing
                elseif eta >= floor(dt/2) && eta <= dt+floor(dt/2)
                    
                    idx_left      = (1:dt:1+n_bar*dt)';
                    idx_right     = (N_r-n_bar*dt:dt:N_r)';
                    idx_dt        = unique(vertcat(idx_left,idx_right));
                    this_route_dt = this_route(idx_dt,2:end);
                    
                    % add an auxiliary point in between eta1 and eta2
                else
                    
                    idx_left            = (1:dt:1+n_bar*dt)';
                    idx_right           = (N_r-n_bar*dt:dt:N_r)';
                    this_route_dt_left  = this_route(idx_left,:);
                    this_route_dt_right = this_route(idx_right,:);
                    info_eta1           = this_route(eta1,:);
                    info_eta2           = this_route(eta2,:);
                    dummy               = zeros(1,length(this_route(1,:)));
                    dummy(1)            = info_eta1(1);
                    dummy(2)            = info_eta1(2);
                    dummy(3)            = info_eta1(3);
                    [lat,lon]           = gcwaypts(info_eta1(4),info_eta1(5),info_eta2(4),info_eta2(5),3);
                    dummy(4)            = lat(2);
                    dummy(5)            = lon(2);
                    dummy(6)            = info_eta1(6);
                    dummy(7)            = info_eta1(7);
                    dummy(8)            = info_eta1(8);
                    dummy(9)            = info_eta1(9);
                    dummy(10)           = info_eta1(10);
                    dummy(11)           = info_eta1(11);
                    
                    this_route_dt       = vertcat(this_route_dt_left(:,2:end),dummy(2:end),this_route_dt_right(:,2:end));
                    
                end
                
                AGG_R_dt     = vertcat(AGG_R_dt,horzcat(cont*ones(numel(this_route_dt(:,1)),1),this_route_dt));
                cont         = cont+1;
                
            end
        else
        end
    end
    
    dlmwrite(strcat(pwd,'/NETWORK/INT_AggRoutes_timestep/',num2str(airport_ID(i)),...
        '/',num2str(airport_ID(i)),'.txt'),AGG_R_dt,' ');
    
    
end
