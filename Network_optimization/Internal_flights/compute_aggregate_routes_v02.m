% Written by Alessandro Bombelli, June 10, 2016
% For each cluster, the aggregate route representative of the cluster is
% computed

clc
clear all
close all

origin      = 34;
destination = 33;

% Checking if the folder exists already or not
if isequal(exist(strcat(pwd,'/NETWORK/INT_AggRoutes'),'dir'),7) % 7 = directory
else
    mkdir(strcat(pwd,'/NETWORK/INT_AggRoutes'));
end

% Checking if the folder exists already or not
if isequal(exist(strcat(pwd,'/NETWORK/INT_AggRoutes/',num2str(origin)),'dir'),7) % 7 = directory
else
    mkdir(strcat(pwd,'/NETWORK/INT_AggRoutes/',num2str(origin)));
end

if exist(strcat(pwd,'/NETWORK/INT_Clusters/',num2str(origin),...
        '/',num2str(origin),'_',num2str(destination),'.txt'), 'file')
    % Loading Cluster matrix
    CLUSTERS = load(strcat(pwd,'/NETWORK/INT_Clusters/',num2str(origin),...
        '/',num2str(origin),'_',num2str(destination),'.txt'));
    % Determining number of clusters (number of rows of the matrix)
    N_clust  = numel(CLUSTERS(:,1));
    % Loading all trajectories
    ALL_TRAJ = load(strcat(pwd,'/NETWORK/INT_traj_filt_info/',num2str(origin),...
        '/',num2str(origin),'_',num2str(destination),'_int.txt'));
    
    AGG_ROUTES = cell(N_clust,1);
    
    for ii=1:N_clust
        
        % Determine number of trajectories for the current cluster
        current_clust_IDs = CLUSTERS(ii,1:numel(find((CLUSTERS(ii,:)~=0)==1)));
        
        % With a singleton trajectory, the trajectory itself is the
        % cluster (we could even decide not to consider the trajectory)
        if numel(current_clust_IDs)==1
            idx              = find(ALL_TRAJ(:,1)==current_clust_IDs);
            agg_route        = ALL_TRAJ(idx,[2:9 16:17]);
            AGG_ROUTES{ii,1} = horzcat(ii*ones(length(agg_route(:,1)),1),agg_route);
        else
            tof = zeros(numel(current_clust_IDs),1);
            for jj=1:numel(current_clust_IDs)
                idx     = find(ALL_TRAJ(:,1)==current_clust_IDs(jj));
                tof(jj) = numel(idx);
            end
            t_low = quantile(tof,0.25);
            t_ar  = round(quantile(tof,0.50));
            t_up  = quantile(tof,0.75);
            idx1  = find(tof<t_low);
            idx2  = find(tof>t_up);
            % Find all trajectories that are compatible with the time of
            % flight
            idx_good = setdiff(1:numel(current_clust_IDs),[idx1;idx2]);
            % At least a trajectory
            if ~isempty(idx_good)
                
                ar_data   = zeros(t_ar,10,numel(idx_good));
                agg_route = zeros(t_ar,10);
                
                for jj=1:numel(idx_good)
                    idx           = find(ALL_TRAJ(:,1)==current_clust_IDs(idx_good(jj)));
                    current_traj  = ALL_TRAJ(idx,:);
                    % Time of flight is already consistent with aggregate
                    % route
                    if numel(idx)==t_ar
                        ar_data(:,:,jj) = current_traj(:,[2:9 16:17]);
                        % Time of flight is smaller than tof of aggregate
                        % route. We need to oversample
                    elseif numel(idx)<t_ar
                        cont = 0;
                        while cont<abs(numel(idx)-t_ar)
                            delta_alt     = abs((current_traj(2:end,8)-current_traj(1:end-1,8))./current_traj(1:end-1,8))*100;
                            idx2          = find(delta_alt<10);
                            dummy         = idx2(2:end)-idx2(1:end-1);
                            dummy2        = dummy(1:end-2)+dummy(2:end-1)+dummy(3:end);
                            idx3          = find(dummy2==3);
                            idx_in        = idx2(idx3(1));
                            idx_fin       = idx2(idx3(end)+3);
                            lat_lon_pairs = current_traj(idx_in:idx_fin,4:5);
                            Hav_dist      = zeros(numel(lat_lon_pairs(:,1))-1,1);
                            for k=1:numel(lat_lon_pairs(:,1))-1
                                Hav_dist(k) = Haversine(lat_lon_pairs(k,:),lat_lon_pairs(k+1,:),1);
                            end
                            B=Hav_dist(1:end-2)+Hav_dist(2:end-1)+Hav_dist(3:end);
                            [~,idx_all]       = sort(B,'descend');
                            idx_addp          = idx_all(1)+idx_in;
                            [add_lat,add_lon] = gcwaypts(current_traj(idx_addp-1,4),current_traj(idx_addp-1,5),...
                                                    current_traj(idx_addp+2,4),current_traj(idx_addp+2,5),5);
                            add_p_latlon      = [add_lat(2:4) add_lon(2:4)];                   
                            addp_data         = current_traj(idx_addp-1:idx_addp+2,:);
                            new_block         = zeros(3,numel(current_traj(1,:)));
                            for k=1:numel(new_block(:,1))
                                new_block(k,1:3)   = addp_data(k,1:3);
                                new_block(k,4:5)   = add_p_latlon(k,:);
                                new_block(k,6)     = mean(addp_data(k:k+1,6));
                                new_block(k,7)     = mean(addp_data(k:k+1,7));
                                new_block(k,8)     = mean(addp_data(k:k+1,8));
                                new_block(k,9)     = mode(addp_data(k:k+1,8));
                                new_block(k,10:15) = addp_data(k,10:15);
                                new_block(k,16)    = mode((addp_data(k:k+1,16)));
                                new_block(k,17)    = mode((addp_data(k:k+1,17)));
                            end
                            
                            new_current_traj = [current_traj(1:idx_addp-1,:);new_block;current_traj(idx_addp+2:end,:)];
                            current_traj     = new_current_traj;
                            cont= cont+1;
                            
                        end
                        ar_data(:,:,jj) = current_traj(:,[2:9 16:17]);
                        % Time of flight is greater than tof of aggregate
                        % route. We need to undersample
                    else
                        cont = 0;
                        while cont<abs(numel(idx)-t_ar)
                            delta_alt     = abs((current_traj(2:end,8)-current_traj(1:end-1,8))./current_traj(1:end-1,8))*100;
                            idx2          = find(delta_alt<10);
                            dummy         = idx2(2:end)-idx2(1:end-1);
                            dummy2        = dummy(1:end-2)+dummy(2:end-1)+dummy(3:end);
                            idx3          = find(dummy2==3);
                            idx_in        = idx2(idx3(1));
                            idx_fin       = idx2(idx3(end)+3);
                            lat_lon_pairs = current_traj(idx_in:idx_fin,4:5);
                            Hav_dist      = zeros(numel(lat_lon_pairs(:,1))-1,1);
                            for k=1:numel(lat_lon_pairs(:,1))-1
                                Hav_dist(k) = Haversine(lat_lon_pairs(k,:),lat_lon_pairs(k+1,:),1);
                            end
                            B=Hav_dist(1:end-2)+Hav_dist(2:end-1)+Hav_dist(3:end);
                            [~,idx_all]       = sort(B,'ascend');
                            idx_addp          = idx_all(1)+idx_in;
                            [add_lat,add_lon] = gcwaypts(current_traj(idx_addp-1,4),current_traj(idx_addp-1,5),...
                                current_traj(idx_addp+2,4),current_traj(idx_addp+2,5),3);
                            add_p_latlon      = [add_lat(2) add_lon(2)];
                            addp_data         = current_traj(idx_addp+1:idx_addp+2,:);
                            new_block         = zeros(1,numel(current_traj(1,:)));
                            for k=1:numel(new_block(:,1))
                                new_block(k,1:3)   = addp_data(k,1:3);
                                new_block(k,4:5)   = add_p_latlon(k,:);
                                new_block(k,6)     = mean(addp_data(k:k+1,6));
                                new_block(k,7)     = mean(addp_data(k:k+1,7));
                                new_block(k,8)     = mean(addp_data(k:k+1,8));
                                new_block(k,9)     = mode(addp_data(k:k+1,8));
                                new_block(k,10:15) = addp_data(k,10:15);
                                new_block(k,16)    = mode((addp_data(k:k+1,16)));
                                new_block(k,17)    = mode((addp_data(k:k+1,17)));
                            end
                            new_current_traj = [current_traj(1:idx_addp-1,:);new_block;current_traj(idx_addp+2:end,:)];
                            current_traj     = new_current_traj;
                            cont= cont+1;
                        end
                        ar_data(:,:,jj) = current_traj(:,[2:9 16:17]);
                    end

                    
                end
                agg_route(:,1)   = ar_data(:,1,1);          % ID of origin airport
                agg_route(:,2)   = ar_data(:,2,1);          % ID of destination airport
                agg_route(:,3)   = mean(ar_data(:,3,:),3);  % latitude
                agg_route(:,4)   = mean(ar_data(:,4,:),3);  % longitude
                agg_route(:,5)   = mean(ar_data(:,5,:),3);  % ground speed
                agg_route(:,6)   = mean(ar_data(:,6,:),3);  % heading
                agg_route(:,7)   = mean(ar_data(:,7,:),3);  % altitude
                agg_route(:,8)   = mode(ar_data(:,8,:),3);  % flight level
                agg_route(:,9)   = mode(ar_data(:,9,:),3);  % Center
                agg_route(:,10)  = mode(ar_data(:,10,:),3); % sector
                AGG_ROUTES{ii,1} = horzcat(ii*ones(length(agg_route(:,1)),1),agg_route);
                
            % No trajectories
            else
            end
        end

    end
    
else
end

AGG_ROUTES_block = [];

for i=1:length(AGG_ROUTES)   
    agg_route = AGG_ROUTES{i,1};
    AGG_ROUTES_block = vertcat(AGG_ROUTES_block,agg_route);
end

dlmwrite(strcat(pwd,'/NETWORK/INT_AggRoutes/',num2str(origin),...
        '/',num2str(origin),'_',num2str(destination),'.txt'),AGG_ROUTES_block,' ');




figure()
hold on
for i=1:length(AGG_ROUTES)   
    agg_route = AGG_ROUTES{i,1};
    plot(agg_route(:,5),agg_route(:,4),'Color','b','Marker','o','Linewidth',2)
end
grid on
axis equal

