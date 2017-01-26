% Written by Alessandro Bombelli, June 10, 2016
% For each cluster, the aggregate route representative of the cluster is
% computed

clc
clear all
close all

origin      = 33;
destination = 34;

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
    
    for ii=1:N_clust
        
        % Determine number of trajectories for the current cluster
        current_clust_IDs = CLUSTERS(ii,1:numel(find((CLUSTERS(ii,:)~=0)==1)));
        
        % With a singleton trajectory, the trajectory itself is the
        % cluster (we could even decide not to consider the trajectory)
        if numel(current_clust_IDs)==1
            idx       = find(ALL_TRAJ(:,1)==current_clust_IDs);
            agg_route = ALL_TRAJ(idx,[2:9 16:17]);
            
        else
            tof = zeros(numel(current_clust_IDs),1);
            for jj=1:numel(current_clust_IDs)
                idx     = find(ALL_TRAJ(:,1)==current_clust_IDs(jj));
                tof(jj) = numel(idx);
            end
            t_low = quantile(tof,0.25);
            t_ar  = quantile(tof,0.50);
            t_up  = quantile(tof,0.75);
            idx1  = find(tof<t_low);
            idx2  = find(tof>t_up);
            % Find all trajectories that are compatible with the time of
            % flight
            idx_good = setdiff(1:numel(current_clust_IDs),[idx1;idx2]);
            % At least a trajectory
            if ~isempty(idx_good)
                for jj=1:numel(idx_good)
                    idx           = find(ALL_TRAJ(:,1)==current_clust_IDs(idx_good(jj)));
                    current_traj  = ALL_TRAJ(idx,:);
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
                    numel(idx)
                    t_ar
                    B=Hav_dist(1:end-2)+Hav_dist(2:end-1)+Hav_dist(3:end)
                    % Time of flight is already consistent with aggregate
                    % route
                    if numel(idx)==t_ar
                    % Time of flight is smaller than tof of aggregate
                    % route. We need to oversample
                    elseif numel(idx)<t_ar
                        add_points  = zeros(abs(t_ar-numel(idx)),8);
                        [~,idx_all] = sort(B,'descend');
                        idx_addp    = idx_all(1:abs(t_ar-numel(idx)))
                        idx_addp+idx_in
                    % Time of flight is greater than tof of aggregate
                    % route. We need to undersample
                    else
                        add_points = zeros(abs(t_ar-numel(idx)),8);
                        [~,idx_all] = sort(B,'ascend');
                        idx_addp    = idx_all(1:abs(t_ar-numel(idx)))
                        idx_addp+idx_in
                    end
                    disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
                end
                
                % No trajectories
            else
            end
        end
        
    disp('%%%%%%%%%%%%%%%%%%%%%%')    
    disp('End of current cluster')
    disp('%%%%%%%%%%%%%%%%%%%%%%')
    end
    
else
end