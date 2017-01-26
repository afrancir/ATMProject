clc
clear all
close all

% Checking if the folder exists already or not
%if isequal(exist(strcat(pwd,'/NETWORK/INT_Clusters'),'dir'),7) % 7 = directory
%else
%    mkdir(strcat(pwd,'/NETWORK/INT_Clusters'));
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loading the Centers that form the planning domain %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid     = fopen(strcat(pwd,'/AUX_DATA/PLANNING_DOMAIN.txt'), 'rt');
CENTERS = textscan(fid,'%s');
fclose(fid);

AIRPORTS    = load(strcat(pwd,'/NETWORK/AIRPORTS/AIRPORTS_coord.txt'));
airport_ID  = AIRPORTS(:,1);

%for i=1:numel(airport_ID)
for i=21
    
    % Checking if the folder exists already or not
    %if isequal(exist(strcat(pwd,'/NETWORK/INT_Clusters/',num2str(airport_ID(i))),'dir'),7) % 7 = directory
    %else
    %    mkdir(strcat(pwd,'/NETWORK/INT_Clusters/',num2str(airport_ID(i))));
    %end
    other_airports = setdiff(airport_ID,airport_ID(i));
    for j=1:numel(other_airports)
        disp(['Trajectories from airport ',num2str(airport_ID(i)),' to airport ',num2str(other_airports(j))]);
        
        latlon_origin_airp = AIRPORTS(airport_ID==airport_ID(i),2:3);
        latlon_dest_airp   = AIRPORTS(airport_ID==other_airports(j),2:3);
        dist_airports      = Haversine(latlon_origin_airp,latlon_dest_airp,1)*6388; % [km]
        
        
        % Check that this OD pair has a valid Frechet Distance Matrix
        if exist(strcat(pwd,'/NETWORK/INT_FDMatrix/',num2str(airport_ID(i)),...
                '/',num2str(airport_ID(i)),'_',num2str(other_airports(j)),'.txt'), 'file')
            
            % Loading all trajectories
            ALL_TRAJ     = load(strcat(pwd,'/NETWORK/INT_traj_filt_info/',num2str(airport_ID(i)),'/',num2str(airport_ID(i)),'_',num2str(other_airports(j)),'_int.txt'));
            flights      = unique(ALL_TRAJ(:,1)); % number of trajectories
            original_IDs = unique(ALL_TRAJ(:,1));
            
            % Loading Frechet Distance Matrix
            FM = load(strcat(pwd,'/NETWORK/INT_FDMatrix/',num2str(airport_ID(i)),...
                '/',num2str(airport_ID(i)),'_',num2str(other_airports(j)),'.txt'));
            
            FM_Almeida = FM;
            
            %%%%%%%%%%%%%%%%%%%%%%%%
            %%% Outliers removal %%%
            %%%%%%%%%%%%%%%%%%%%%%%%
            
            n_deleted = 1;
            outliers  = [];
            c_lim     = 0.03;
            
            while n_deleted>0;
                
                full_FM        = FM_Almeida+FM_Almeida';
                full_FM_sorted = sort(full_FM,2);
                d              = full_FM_sorted(:,2);
                N              = numel(flights);
                avg_d          = sum(d)/numel(d);
                R              = 4*avg_d;
                [~,col_ids]    = find(full_FM<R);
                c = zeros(N,1);
                for k=1:N;
                    c(k) = numel(find(col_ids==k))-1;
                end
                avg_c = sum(c)/N;                   % average connectivity
                idxs = find(c<c_lim*avg_c);         % routes to delete in this first iteration
                outliers=[outliers; flights(idxs)]; % we save all the outliers
                n_deleted=numel(idxs);
                % We delete the flights with the furthest neighbors
                FM_Almeida(idxs,:)  = [];
                FM_Almeida(:, idxs) = [];
                flights(idxs)       = [];
            end
            
            
            n_deleted = 1;
            while n_deleted>0;
                
                full_FM        = FM_Almeida+FM_Almeida';
                full_FM_sorted = sort(full_FM,2);
                d              = full_FM_sorted(:,2);
                N              = numel(flights);
                avg_d          = sum(d)/numel(d);
                R              = 2*avg_d;
                [~,col_ids]    = find(full_FM<R);
                c = zeros(N,1);
                for k=1:N;
                    c(k) = numel(find(col_ids==k))-1;
                end
                avg_c = sum(c)/N;                   % average connectivity
                idxs = find(c<c_lim*avg_c);         % routes to delete in this first iteration
                outliers=[outliers; flights(idxs)]; % we save all the outliers
                n_deleted=numel(idxs);
                % We delete the flights with the furthest neighbors
                FM_Almeida(idxs,:) = [];
                FM_Almeida(:,idxs) = [];
                flights(idxs)      = [];
            end
            
            [IDX,isnoise] = DBSCAN_v01(FM,0.1*dist_airports,2);
            
            if numel(outliers) ~= 0
                outliers_Almeida = outliers;
            else
                outliers_Almeida = NaN;
            end
            
            if ~isempty(find(isnoise,1))
                outliers_DBSCAN = original_IDs(isnoise);
            else
                outliers_DBSCAN = NaN;
            end
            
            %disp(['Outliers from airport ',num2str(airport_ID(i)),' to airport ',num2str(other_airports(j)),':'])
            %outliers_Almeida
            %outliers_DBSCAN
            
            outliers                  = union(outliers_Almeida,outliers_DBSCAN);
            outliers(isnan(outliers)) = [];
            
            % flight_ID_wo_outliers: IDs of the flight tracks we want to
            % keep for the clustering process
            % idx_flight_ID_wo_outliers: indexes of the flights tracks we
            % want to keep, i.e., position of these tracks in the Frechet
            % distance matrix
            [flight_ID_wo_outliers,idx_flight_ID_wo_outliers] = setdiff(original_IDs,outliers);
            % idx_outliers: indexes of the outliers, i.e.,  position of
            % the outliers in the Frechet distance matrix
            idx_outliers                                      = setdiff(1:numel(original_IDs),idx_flight_ID_wo_outliers);
            
            disp(['Number of trajectories: ',num2str(numel(flight_ID_wo_outliers))])
            
            
            if numel(flight_ID_wo_outliers)>=2
                
                rgb_airports = 0.5*ones(3,1);
                rgb_all_traj = 0.7*ones(3,1);
                rgb_outliers = 0.3*ones(3,1);
                
                figure()
                hold on
                
                % Plot outliers, if any
                if ~isnan(outliers)
                    for ii=1:numel(outliers)
                        idx = find(ALL_TRAJ(:,1)==outliers(ii));
                        plot(ALL_TRAJ(idx,5),ALL_TRAJ(idx,4),'Color',rgb_outliers,'Linestyle','none','Linewidth',1,'Marker','*','Markersize',2.5)
                    end
                else
                end
                
                % Plot "good flights", if any
                if ~isnan(flight_ID_wo_outliers)
                    for ii=1:numel(flight_ID_wo_outliers)
                        idx = find(ALL_TRAJ(:,1)==flight_ID_wo_outliers(ii));
                        plot(ALL_TRAJ(idx,5),ALL_TRAJ(idx,4),'Color',rgb_all_traj,'Linestyle','none','Linewidth',1,'Marker','*','Markersize',2.5)
                    end
                else
                end
                
                plot(AIRPORTS(:,3),AIRPORTS(:,2),'Color',rgb_airports,'Linewidth',2,'Linestyle','none','Marker','s','Markersize',10)
                for ii=1:length(CENTERS{1,1})
                    current_Center_latlon = load(strcat(pwd,'/NETWORK/Center_boundaries/',char(CENTERS{1,1}(ii)),'.txt'));
                    plot(current_Center_latlon(:,2),current_Center_latlon(:,1),'Color','k','Linewidth',2,'Marker','none')
                end
                grid on
                axis equal
                xlabel('Longitude [deg]','Fontname','Avantgarde','Fontsize',14)
                ylabel('Latiitude [deg]','Fontname','Avantgarde','Fontsize',14)
                
                % Eliminate rows and columns in the Frechet distance matrix
                % that refer to outliers
                
                FM(idx_outliers,:) = [];
                FM(:,idx_outliers) = [];
                
                N_traj      = length(FM(:,1));
                dist_mat_km = FM+FM';
                dist_vec_km = squareform(dist_mat_km);
                tree        = linkage(dist_vec_km,'average');
                
                LB          = min(30,0.05*dist_airports); % [km]
                UB          = max(90,0.15*dist_airports); % [km]
                idx_LB      = find(tree(:,3)>=LB);
                idx_UB      = find(tree(:,3)<=UB);
                if numel(idx_LB)>=1
                    idx_LB = [idx_LB(1)-1;idx_LB];
                else
                end
                idx         = intersect(idx_LB,idx_UB);
                N_c         = sort(N_traj-idx);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 h = dendrogram(tree,30,'ColorThreshold',0.4*max(tree(:,3)));
%                 %set(h,'Position',pos)
%                 dummy_gray = zeros(numel(h),3);
%                 for jj=1:numel(h)
%                     dummy_gray(jj,:) = (0.1+0.3*rand(1))*ones(1,3);
%                 end
%                 cm=num2cell(dummy_gray,2);
%                 set(h,...
%                     'marker','s',...
%                     'markersize',6,...
%                     {'markerfacecolor'},num2cell(ones(numel(h),3),2),...
%                     {'markeredgecolor'},cm,...
%                     'linewidth',1.5,...
%                     {'color'},cm);
%                 xlabel('Flight track ID','Fontname','Avantgarde','Fontsize',14)
%                 ylabel('Distance [km]','Fontname','Avantgarde','Fontsize',14)
%                 set(h,'LineWidth',2)
%                 xLimits = get(gca,'XLim');
%                 yLimits = get(gca,'YLim');
%                 hold on
%                 plot(xLimits,[LB LB],'Color',0.3*ones(1,3),'Linestyle','--','Linewidth',2)
%                 plot(xLimits,[UB UB],'Color',0.5*ones(1,3),'Linestyle','--','Linewidth',2)
%                 text(5,50,'Lower boundary','Fontname','Avantgarde','Fontsize',14)
%                 %text(5,600,'Upper boundary','Fontname','Avantgarde','Fontsize',14)
%                 grid on
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if numel(N_c) == 1
                    NC_opt = N_c;
                    T_opt  = cluster(tree,'maxclust',NC_opt);
                elseif numel(N_c) == 0
                    NC_opt = 1;
                    T_opt  = cluster(tree,'maxclust',NC_opt);
                else
                    if N_c(1) == 1
                        N_c = [N_c(2:end);N_c(end)+1];
                        % At most, N_c(end) = N_traj-1, thus N_c(end)+1
                        % will be N_traj in the worst case scenario (it
                        % corresponds to the case where each singleton
                        % forms its own cluster)
                    else
                        N_c = [N_c;N_c(end)+1];
                    end
                    
                    AS   = zeros(numel(N_c),1);
                    DB   = zeros(numel(N_c),1);
                    DI   = zeros(numel(N_c),1);
                    
                    flag_indices = 0;
                    
                    for k=1:numel(N_c)
                        
                        N_clust    = N_c(k);
                        T          = cluster(tree,'maxclust',N_clust);
                        
                        s      = silhouette([],T,dist_vec_km);
                        AS(k)  = mean(s);
                        DB(k)  = DaviesBouldin(dist_mat_km,N_clust,T,flag_indices);
                        DI(k)  = Dunn(dist_mat_km,N_clust,T,flag_indices);
                        
                    end
                    
                    if N_c(1) == 2
                        N_c_aug = [1;N_c];
                        AS_aug  = [0;AS];
                        DB_aug  = [1;DB];
                        DI_aug  = [1;DI];
                    else
                        N_c_in  = N_c(1)-1;
                        T_in    = cluster(tree,'maxclust',N_c_in);
                        s_in    = silhouette([],T_in,dist_vec_km);
                        AS_in   = mean(s_in);
                        DB_in   = DaviesBouldin(dist_mat_km,N_c_in,T_in,flag_indices);
                        DI_in   = Dunn(dist_mat_km,N_c_in,T_in,flag_indices);
                        N_c_aug = [N_c_in;N_c];
                        AS_aug  = [AS_in;AS];
                        DB_aug  = [DB_in;DB];
                        DI_aug  = [DI_in;DI];
                    end
                    
                    local_maxmin    = zeros(3,1);
                    NC_optimal      = zeros(3,1);
                    
                    [val_AS,idx_AS] = findpeaks(AS_aug);
                    if isempty(val_AS)
                    else
                        local_maxmin(1,1) = 1;
                        if numel(val_AS(:,1))>1
                            [~,dummy]       = max(val_AS(:,1));
                            NC_optimal(1,1) = N_c_aug(idx_AS(dummy));
                        else
                            NC_optimal(1,1) = N_c_aug(idx_AS);
                        end
                    end
                    [val_DB,idx_DB] = findpeaks(-DB_aug);
                    if isempty(val_DB)
                    else
                        local_maxmin(2,1) = 1;
                        if numel(val_DB(:,1))>1
                            [~,dummy]       = max(val_DB(:,1));
                            NC_optimal(2,1) = N_c_aug(idx_DB(dummy));
                        else
                            NC_optimal(2,1) = N_c_aug(idx_DB);
                        end
                    end
                    [val_DI,idx_DI] = findpeaks(DI_aug);
                    if isempty(val_DI)
                    else
                        local_maxmin(3,1) = 1;
                        if numel(val_DI(:,1))>1
                            [~,dummy]       = max(val_DI(:,1));
                            NC_optimal(3,1) = N_c_aug(idx_DI(dummy));
                        else
                            NC_optimal(3,1) = N_c_aug(idx_DI);
                        end
                    end
                    
                    idx_with_maxmin = find(local_maxmin);
                    
                    if isempty(find(local_maxmin==1,1))
                        NC_opt = 1;
                    else
                        NC_opt = min(NC_optimal(idx_with_maxmin,1));
                        
                        rgb1   = 0.1*ones(1,3);
                        rgb2   = 0.2*ones(1,3);
                        rgb3   = 0.3*ones(1,3);
                        
                        idx = find(N_c_aug==NC_opt);
                        
%                         hh = figure();
%                         subplot(2,1,1)
%                         hold on
%                         plot(N_c_aug,AS_aug,'Color',rgb1,'Linestyle','-','Linewidth',2,'Marker','s','Markersize',6)
%                         plot(N_c_aug,DB_aug,'Color',rgb2,'Linestyle','-','Linewidth',2,'Marker','*','Markersize',6)
%                         yLimits = get(gca,'YLim');
%                         plot([N_c_aug(2)-0.1 N_c_aug(2)-0.1],yLimits,'Color','k','Linestyle','--','Linewidth',2)
%                         plot([N_c_aug(end-1)+0.1 N_c_aug(end-1)+0.1],yLimits,'Color','k','Linestyle','--','Linewidth',2)
%                         plot(N_c_aug(idx),AS_aug(idx),'Color','k','Linestyle','-','Linewidth',4,'Marker','+','Markersize',12)
%                         plot(N_c_aug(idx),DB_aug(idx),'Color','k','Linestyle','-','Linewidth',4,'Marker','+','Markersize',12)
%                         h_legend1 = legend('Average Silhouette','Davies-Bouldin Index','Location','Best');
%                         set(h_legend1,'Fontname','Avantgarde','Fontsize',10)
%                         set(gca,'Ytick',[0 1])
%                         ylabel('Performance index','Fontname','Avantgarde','Fontsize',14)
%                         grid on
%                         subplot(2,1,2)
%                         hold on
%                         plot(N_c_aug,DI_aug,'Color',rgb3,'Linestyle','-','Linewidth',2,'Marker','o','Markersize',6)
%                         yLimits = get(gca,'YLim');
%                         plot([N_c_aug(2)-0.1 N_c_aug(2)-0.1],yLimits,'Color','k','Linestyle','--','Linewidth',2)
%                         plot([N_c_aug(end-1)+0.1 N_c_aug(end-1)+0.1],yLimits,'Color','k','Linestyle','--','Linewidth',2)
%                         plot(N_c_aug(idx),DI_aug(idx),'Color','k','Linestyle','-','Linewidth',4,'Marker','+','Markersize',12)
%                         h_legend2 = legend('Dunn Index','Location','Best');
%                         set(h_legend2,'Fontname','Avantgarde','Fontsize',10)
%                         xlabel('Number of clusters','Fontname','Avantgarde','Fontsize',14)
%                         ylabel('Performance index','Fontname','Avantgarde','Fontsize',14)
%                         grid on
%                         
                        T_opt = cluster(tree,'maxclust',NC_opt);
                        
                        c = zeros(NC_opt,3);
                        
                        for jj=1:NC_opt
                            c(jj,:) = (0.1+0.6*rand(1))*ones(1,3);
                        end
                        
%                         figure()
%                         hold on
%                         for ii=1:NC_opt
%                             position_tracks = find(T_opt==ii);
%                             tracks_ID       = flight_ID_wo_outliers(position_tracks);
%                             for jj=1:numel(tracks_ID)
%                                 this_track_ID  = tracks_ID(jj);
%                                 this_track_idx = find(ALL_TRAJ(:,1)==this_track_ID);
%                                 plot(ALL_TRAJ(this_track_idx,5),ALL_TRAJ(this_track_idx,4),'Color',c(ii,:),'Linestyle','none','Linewidth',1,'Marker','*','Markersize',2.5)
%                             end
%                         end
%                         plot(AIRPORTS(:,3),AIRPORTS(:,2),'Color',rgb_airports,'Linewidth',2,'Linestyle','none','Marker','s','Markersize',10)
%                         for ii=1:length(CENTERS{1,1})
%                             current_Center_latlon = load(strcat(pwd,'/NETWORK/Center_boundaries/',char(CENTERS{1,1}(ii)),'.txt'));
%                             plot(current_Center_latlon(:,2),current_Center_latlon(:,1),'Color','k','Linewidth',2,'Marker','none')
%                         end
%                         grid on
%                         axis equal
%                         xlabel('Longitude [deg]','Fontname','Avantgarde','Fontsize',14)
%                         ylabel('Latiitude [deg]','Fontname','Avantgarde','Fontsize',14)
                        
                    end
                end
                
                CLUSTERS = zeros(NC_opt,N_traj);
                for ii=1:NC_opt
                    idx = find(T_opt==ii);
                    if ii==1
                        max_traj = numel(idx);
                    else
                        if numel(idx)>max_traj
                            max_traj = numel(idx);
                        else
                        end
                    end
                    CLUSTERS(ii,1:numel(idx)) = transpose(flights(idx));
                end
                % Getting rid of extra columns
                CLUSTERS = CLUSTERS(:,1:max_traj);
                
                
                
                
                % Storing matrix in the associated folder
                dlmwrite(strcat(pwd,'/NETWORK/INT_Clusters/',num2str(airport_ID(i)),'/',num2str(airport_ID(i)),'_',num2str(other_airports(j)),'.txt'),CLUSTERS,' ');
                disp(['Number of clusters: ',num2str(NC_opt)]);
                disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
                
                
                
                
            else
            end
            
            
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            
            
        else
            disp(['No valid trajectories from ',num2str(airport_ID(i)),' to airport ',num2str(other_airports(j))])
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        end
    end
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('%%% End current airport%%%')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%')
end

            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
