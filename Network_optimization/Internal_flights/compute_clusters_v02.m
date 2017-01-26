clc
clear all
close all

% Checking if the folder exists already or not
if isequal(exist(strcat(pwd,'/NETWORK/INT_Clusters'),'dir'),7) % 7 = directory
else
    mkdir(strcat(pwd,'/NETWORK/INT_Clusters'));
end

% Checking if the folder exists already or not
if isequal(exist(strcat(pwd,'/FIGURES'),'dir'),7) % 7 = directory
else
    mkdir(strcat(pwd,'/FIGURES'));
end

AIRPORTS    = load(strcat(pwd,'/NETWORK/AIRPORTS/AIRPORTS_coord.txt'));
airport_ID = AIRPORTS(:,1);

for i=1:numel(airport_ID)
    % Checking if the folder exists already or not
    if isequal(exist(strcat(pwd,'/NETWORK/INT_Clusters/',num2str(airport_ID(i))),'dir'),7) % 7 = directory
    else
        mkdir(strcat(pwd,'/NETWORK/INT_Clusters/',num2str(airport_ID(i))));
    end
    other_airports = setdiff(airport_ID,airport_ID(i));
    for j=1:numel(other_airports)
        disp(['Trajectories from airport ',num2str(airport_ID(i)),' to airport ',num2str(other_airports(j))]);
        
        % Check that this OD pair has a valid Frechet Distance Matrix
        if exist(strcat(pwd,'/NETWORK/INT_FDMatrix/',num2str(airport_ID(i)),...
                '/',num2str(airport_ID(i)),'_',num2str(other_airports(j)),'.txt'), 'file')
            
            % Loading all trajectories
            ALL_TRAJ = load(strcat(pwd,'/NETWORK/INT_traj_filt_info/',num2str(airport_ID(i)),'/',num2str(airport_ID(i)),'_',num2str(other_airports(j)),'_int.txt'));
            flights = unique(ALL_TRAJ(:,1)); % number of trajectories
            
            % Loading Frechet Distance Matrix
            FM = load(strcat(pwd,'/NETWORK/INT_FDMatrix/',num2str(airport_ID(i)),...
                '/',num2str(airport_ID(i)),'_',num2str(other_airports(j)),'.txt'));
            
            %%%%%%%%%%%%%%%%%%%%%%%%
            %%% Outliers removal %%%
            %%%%%%%%%%%%%%%%%%%%%%%%
            
            n_deleted = 1;
            outliers  = [];
            c_lim     = 0.03;
            
            while n_deleted>0;
                
                full_FM        = FM+FM';
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
                FM(idxs,:)=[];
                FM(:, idxs)=[];
                flights(idxs)=[];
            end
            
            
            n_deleted = 1;
            while n_deleted>0;
                
                full_FM        = FM+FM';
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
                FM(idxs,:)=[];
                FM(:, idxs)=[];
                flights(idxs)=[];
            end
            disp(['Number of outliers removed: ',num2str(numel(outliers))]);
            disp(['Number of trajectories remaining: ',num2str(numel(flights))]);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Selection of the optimal number of clusters %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            R_e         = 6371;  % Earth's radius [km]
            h           = 10;    % altitude [km]
            radius      = R_e+h; % radius of the sphere where we project the points
            
            N_traj = numel(flights);
            
            % If we have only 2 trajectories characterizing this OD pair,
            % we check whether the onlt non-zero component of the Frechet
            % Distance Matrix is bigger or smaller than a specific
            % threshold. In the first case, we consider each trajectory as
            % a single cluster. In the second case, the two trajectories
            % form a unique cluster
            if N_traj==2
                origin_airp_idx = find(AIRPORTS(:,1)==airport_ID(i));
                dest_airp_idx   = find(AIRPORTS(:,1)==other_airports(j));
                dist_airp       = Haversine(AIRPORTS(origin_airp_idx,2:3),AIRPORTS(dest_airp_idx,2:3),1)*radius;
                
                % In this case, we consider the two trajectories
                % sufficiently close to belong to the same cluster
                if FM(1,1)<0.4*dist_airp
                    T        = [1;1];
                    CLUSTERS = flights';
                % In this case, let's split     
                else
                    T        = [1;2];
                    CLUSTERS = flights;
                end
                
                
            else
                origin_airp_idx = find(AIRPORTS(:,1)==airport_ID(i));
                dest_airp_idx   = find(AIRPORTS(:,1)==other_airports(j));
                dist_airp       = Haversine(AIRPORTS(origin_airp_idx,2:3),AIRPORTS(dest_airp_idx,2:3),1)*radius;
                mean_FD         = mean(FM(FM~=0));
                
                % Reshape dist_mat_km into a vector compatible with the matlab's pdist
                % format
                dist_vec_km = [];
                for ii=1:N_traj-1
                    dist_vec_km= horzcat(dist_vec_km,FM(ii,ii+1:end));
                end
                
                max_clusters = 10;
                N_c          = 2:min(N_traj,max_clusters);
                
               
                dist_mat_km = FM;
                
                s_mean_vec    = zeros(numel(N_c),1);
                s_av_mean_vec = zeros(numel(N_c),1);
                DB_vec        = zeros(numel(N_c),1);
                DB_av_vec     = zeros(numel(N_c),1);
                
                for k=1:numel(N_c)
                    
                    N_clust    = N_c(k);
                    dendrogram = linkage(dist_vec_km,'average');
                    T          = cluster(dendrogram,'maxclust',N_clust);
                    
                    s = silhouette([],T,dist_vec_km);
                    
                    s_vec        = zeros(N_clust,1);
                    num_traj_vec = zeros(N_clust,1);
                    for ii=1:N_clust
                        idx             = find(T==ii);
                        num_traj_vec(ii) = numel(idx);
                        s_vec(ii)        = mean(s(idx));
                    end
                    
                    s_mean_vec(k)    = mean(s_vec);
                    s_av_mean_vec(k) = dot(num_traj_vec,s_vec)/numel(T);
                    
                end
                
                % get complete Frechet distance matrix that will be used when computing
                % extra-cluster distance
                FM_complete = FM+FM';
                
                for k=1:numel(N_c)
                    
                    N_clust    = N_c(k);
                    dendrogram = linkage(dist_vec_km,'average');
                    T          = cluster(dendrogram,'maxclust',N_clust);
                    
                    num_traj_vec = zeros(N_clust,1);
                    for ii=1:N_clust
                        idx             = find(T==ii);
                        num_traj_vec(ii) = numel(idx);
                    end
                    
                    intra_cluster_distance = zeros(N_clust,1);
                    extra_cluster_distance = zeros(N_clust,N_clust);
                    R                      = zeros(N_clust,N_clust);
                    D                      = zeros(N_clust,1);
                    
                    % Compute all intra-cluster distances
                    for ii=1:N_clust
                        % find rows/columns associated with the current cluster
                        r_c  = find(T==ii);
                        % if cluster has only one element, set intra-cluster distance to
                        % zero
                        if numel(r_c) == 1
                            intra_cluster_distance(ii) = 0;
                        else
                            % isolate sub-matrix associated with current cluster
                            FM_local = FM(r_c,r_c);
                            % compute average intra-cluster distance
                            intra_cluster_distance(ii) = mean(FM_local(FM_local~=0));
                        end
                    end
                    
                    % Compute all extra-cluster distances
                    for ii=1:N_clust
                        for jj=1:N_clust
                            if jj==ii
                            else
                                rows     = find(T==ii);
                                cols     = find(T==jj);
                                FM_local = FM_complete(rows,cols);
                                extra_cluster_distance(ii,jj) = mean(FM_local(FM_local~=0));
                            end
                        end
                    end
                    
                    for ii=1:N_clust
                        for jj=1:N_clust
                            if ii==jj
                            else
                                R(ii,jj) = (intra_cluster_distance(ii)+intra_cluster_distance(jj))/...
                                    (extra_cluster_distance(ii,jj));
                            end
                        end
                    end
                    
                    for ii=1:N_clust
                        D(ii) = max(R(ii,1:end ~= ii));
                    end
                    DB_vec(k)    = sum(D)/N_clust;
                    DB_av_vec(k) = dot(num_traj_vec,D)/numel(T);
                end
                
                fake_point_DB = 2;
                fake_point_AS = 0;
                
                s_mean_vec    = vertcat(fake_point_AS,s_mean_vec);
                s_av_mean_vec = vertcat(fake_point_AS,s_av_mean_vec);
                DB_vec        = vertcat(fake_point_DB,DB_vec);
                DB_av_vec     = vertcat(fake_point_DB,DB_av_vec);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Find local maximum for average silhouette %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if numel(s_mean_vec) == 3;
                    [~,idx_AS] = max(s_mean_vec);
                    flag_AS    = 0;
                else
                    
                    [value_AS,idx_AS] = findpeaks(s_mean_vec);
                    % No local maximum was identified. In this case, we
                    % consider a single cluster as the best option
                    if isempty(idx_AS)
                        flag_AS = 1;
                    else
                        flag_AS      = 0;
                        [~,dummy_AS] = max(value_AS);
                        idx_AS       = idx_AS(dummy_AS);
                    end
                    
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Find local maximum for weighted average silhouette %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if numel(s_av_mean_vec) == 3;
                    [~,idx_WAS] = max(s_mean_vec);
                    flag_WAS    = 0;
                else
                    
                    [value_WAS,idx_WAS] = findpeaks(s_av_mean_vec);
                    % No local maximum was identified. In this case, we
                    % consider a single cluster as the best option
                    if isempty(idx_WAS)
                        flag_WAS = 1;
                    else
                        flag_WAS      = 0;
                        [~,dummy_WAS] = max(value_WAS);
                        idx_WAS       = idx_WAS(dummy_WAS);
                    end
                    
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Find local minimum for DB index %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if numel(DB_vec) == 3;
                    [~,idx_DB] = min(DB_vec);
                    flag_DB    = 0;
                else
                    
                    % Since we are looking for a local minimum in this case, we
                    % search for local maxima of -DB_vec
                    [value_DB,idx_DB] = findpeaks(-DB_vec);
                    % No local maximum was identified. In this case, we
                    % consider a single cluster as the best option
                    if isempty(idx_DB)
                        flag_DB = 1;
                    else
                        flag_DB      = 0;
                        [~,dummy_DB] = max(value_DB);
                        idx_DB       = idx_DB(dummy_DB);
                    end
                    
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Find local minimum for weighted DB index %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if numel(DB_av_vec) == 3;
                    [~,idx_WDB] = min(DB_av_vec);
                    flag_DB    = 0;
                else
                    % Since we are looking for a local minimum in this case, we
                    % search for local maxima of -WDB_vec
                    [value_WDB,idx_WDB] = findpeaks(-DB_av_vec);
                    % No local maximum was identified. In this case, we
                    % consider a single cluster as the best option
                    if isempty(idx_WDB)
                        flag_WDB = 1;
                    else
                        flag_WDB      = 0;
                        [~,dummy_WDB] = max(value_WDB);
                        idx_WDB       = idx_WDB(dummy_WDB);
                    end
                    
                end
                    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Now compare the four indexes and determine what is %%%
                %%% the optimal number of clusters                     %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                flag_vector = [flag_DB;flag_WDB;flag_AS;flag_WAS];
                
                % If at least 2 indexes do not show any local
                % maximum/minimum, the set of trajectories form a single
                % cluster
                if numel(find(flag_vector))>=2
                    
                    N_clusters_selected = 1;
                    T                   = ones(N_traj,1);
                    % Otherwise, discard the index that is not showing a local
                    % maximum/minimum (if present), and select the number of
                    % clusters checking the best overall performance
                else
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%% All indexes show a minimum/maximum %%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if numel(find(flag_vector))==0
                        
                        % if the four indexes match, any of them can be chosen as the optimal
                        % number of clusters
                        idx_vec = [idx_DB,idx_WDB,idx_AS,idx_WAS];
                        
                        if (numel(unique(idx_vec))==1)
                            idx_min = idx_DB;
                            % otherwise, choose the one that minimizes the relative error of the other
                            % indexes
                        else
                            index_err_matrix      = zeros(4,4);
                            index_err_matrix(1,2) = abs((DB_vec(idx_WDB)-DB_vec(idx_DB))/DB_vec(idx_DB))*100;
                            index_err_matrix(1,3) = abs((DB_vec(idx_AS)-DB_vec(idx_DB))/DB_vec(idx_DB))*100;
                            index_err_matrix(1,4) = abs((DB_vec(idx_WAS)-DB_vec(idx_DB))/DB_vec(idx_DB))*100;
                            index_err_matrix(2,1) = abs((DB_av_vec(idx_DB)-DB_av_vec(idx_WDB))/DB_av_vec(idx_WDB))*100;
                            index_err_matrix(2,3) = abs((DB_av_vec(idx_AS)-DB_av_vec(idx_WDB))/DB_av_vec(idx_WDB))*100;
                            index_err_matrix(2,4) = abs((DB_av_vec(idx_WAS)-DB_av_vec(idx_WDB))/DB_av_vec(idx_WDB))*100;
                            index_err_matrix(3,1) = abs((s_mean_vec(idx_DB)-s_mean_vec(idx_AS))/s_mean_vec(idx_AS))*100;
                            index_err_matrix(3,2) = abs((s_mean_vec(idx_WDB)-s_mean_vec(idx_AS))/s_mean_vec(idx_AS))*100;
                            index_err_matrix(3,4) = abs((s_mean_vec(idx_WAS)-s_mean_vec(idx_AS))/s_mean_vec(idx_AS))*100;
                            index_err_matrix(4,1) = abs((s_av_mean_vec(idx_DB)-s_av_mean_vec(idx_WAS))/s_av_mean_vec(idx_WAS))*100;
                            index_err_matrix(4,2) = abs((s_av_mean_vec(idx_WDB)-s_av_mean_vec(idx_WAS))/s_av_mean_vec(idx_WAS))*100;
                            index_err_matrix(4,3) = abs((s_av_mean_vec(idx_AS)-s_av_mean_vec(idx_WAS))/s_av_mean_vec(idx_WAS))*100;
                            
                            err_vec               = [mean(index_err_matrix(1:end ~=1,1)),mean(index_err_matrix(1:end ~=2,2)),...
                                mean(index_err_matrix(1:end ~=3,3)),mean(index_err_matrix(1:end ~=4,4))];
                            [~,idx_err_vec] = min(err_vec);
                            idx_min         = idx_vec(idx_err_vec);
                        end
                      
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
                    %%% One index is discarded %%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    else
                        % DB index is discarded
                        if find(flag_vector)==1
                            
                            idx_vec               = [idx_WDB,idx_AS,idx_WAS];
                            
                            index_err_matrix      = zeros(3,3);
                            index_err_matrix(1,2) = abs((DB_av_vec(idx_WDB)-DB_av_vec(idx_AS))/DB_av_vec(idx_WDB))*100;
                            index_err_matrix(1,3) = abs((DB_av_vec(idx_WDB)-DB_av_vec(idx_WAS))/DB_av_vec(idx_WDB))*100;
                            index_err_matrix(2,1) = abs((s_mean_vec(idx_AS)-s_mean_vec(idx_WDB))/s_mean_vec(idx_AS))*100;
                            index_err_matrix(2,3) = abs((s_mean_vec(idx_AS)-s_mean_vec(idx_WAS))/s_mean_vec(idx_AS))*100;
                            index_err_matrix(3,1) = abs((s_av_mean_vec(idx_WAS)-s_av_mean_vec(idx_WDB))/s_av_mean_vec(idx_WAS))*100;
                            index_err_matrix(3,2) = abs((s_av_mean_vec(idx_WAS)-s_av_mean_vec(idx_AS))/s_av_mean_vec(idx_WAS))*100;
                            
                            err_vec               = [mean(index_err_matrix(1:end ~=1,1)),mean(index_err_matrix(1:end ~=2,2)),...
                                mean(index_err_matrix(1:end ~=3,3))];
                            [~,idx_err_vec] = min(err_vec);
                            idx_min         = idx_vec(idx_err_vec);
                        % weighted DB index is discarded
                        elseif find(flag_vector)==2
                            
                            idx_vec               = [idx_DB,idx_AS,idx_WAS];
                            
                            index_err_matrix      = zeros(3,3);
                            index_err_matrix(1,2) = abs((DB_vec(idx_DB)-DB_vec(idx_AS))/DB_vec(idx_DB))*100;
                            index_err_matrix(1,3) = abs((DB_vec(idx_DB)-DB_vec(idx_WAS))/DB_vec(idx_DB))*100;
                            index_err_matrix(2,1) = abs((s_mean_vec(idx_AS)-s_mean_vec(idx_DB))/s_mean_vec(idx_AS))*100;
                            index_err_matrix(2,3) = abs((s_mean_vec(idx_AS)-s_mean_vec(idx_WAS))/s_mean_vec(idx_AS))*100;
                            index_err_matrix(3,1) = abs((s_av_mean_vec(idx_WAS)-s_av_mean_vec(idx_DB))/s_av_mean_vec(idx_WAS))*100;
                            index_err_matrix(3,2) = abs((s_av_mean_vec(idx_WAS)-s_av_mean_vec(idx_AS))/s_av_mean_vec(idx_WAS))*100;
                            
                            err_vec               = [mean(index_err_matrix(1:end ~=1,1)),mean(index_err_matrix(1:end ~=2,2)),...
                                mean(index_err_matrix(1:end ~=3,3))];
                            [~,idx_err_vec] = min(err_vec);
                            idx_min         = idx_vec(idx_err_vec);
                        % average silhouette is discarded
                        elseif find(flag_vector)==3
                            
                            idx_vec               = [idx_DB,idx_WDB,idx_WAS];
                            
                            index_err_matrix      = zeros(3,3);
                            index_err_matrix(1,2) = abs((DB_vec(idx_DB)-DB_vec(idx_WDB))/DB_vec(idx_DB))*100;
                            index_err_matrix(1,3) = abs((DB_vec(idx_DB)-DB_vec(idx_WAS))/DB_vec(idx_DB))*100;
                            index_err_matrix(2,1) = abs((DB_av_vec(idx_WDB)-DB_av_vec(idx_DB))/DB_av_vec(idx_WDB))*100;
                            index_err_matrix(2,3) = abs((DB_av_vec(idx_WDB)-DB_av_vec(idx_WAS))/DB_av_vec(idx_WDB))*100;
                            index_err_matrix(3,1) = abs((s_av_mean_vec(idx_WAS)-s_av_mean_vec(idx_DB))/s_av_mean_vec(idx_WAS))*100;
                            index_err_matrix(3,2) = abs((s_av_mean_vec(idx_WAS)-s_av_mean_vec(idx_WDB))/s_av_mean_vec(idx_WAS))*100;
                            
                            err_vec               = [mean(index_err_matrix(1:end ~=1,1)),mean(index_err_matrix(1:end ~=2,2)),...
                                mean(index_err_matrix(1:end ~=3,3))];
                            [~,idx_err_vec] = min(err_vec);
                            idx_min         = idx_vec(idx_err_vec);
                        % weighted average silhouette is discarded
                        else
                            
                            idx_vec               = [idx_DB,idx_WDB,idx_AS];
                            
                            index_err_matrix      = zeros(3,3);
                            index_err_matrix(1,2) = abs((DB_vec(idx_DB)-DB_vec(idx_WDB))/DB_vec(idx_DB))*100;
                            index_err_matrix(1,3) = abs((DB_vec(idx_DB)-DB_vec(idx_AS))/DB_vec(idx_DB))*100;
                            index_err_matrix(2,1) = abs((DB_av_vec(idx_WDB)-DB_av_vec(idx_DB))/DB_av_vec(idx_WDB))*100;
                            index_err_matrix(2,3) = abs((DB_av_vec(idx_WDB)-DB_av_vec(idx_AS))/DB_av_vec(idx_WDB))*100;
                            index_err_matrix(3,1) = abs((s_mean_vec(idx_AS)-s_mean_vec(idx_DB))/s_mean_vec(idx_AS))*100;
                            index_err_matrix(3,2) = abs((s_mean_vec(idx_AS)-s_mean_vec(idx_WDB))/s_mean_vec(idx_AS))*100;
                            
                            err_vec               = [mean(index_err_matrix(1:end ~=1,1)),mean(index_err_matrix(1:end ~=2,2)),...
                                mean(index_err_matrix(1:end ~=3,3))];
                            [~,idx_err_vec] = min(err_vec);
                            idx_min         = idx_vec(idx_err_vec);
                        end
                        
                    end
                    
                    N_clusters_selected = N_c(idx_min);
                    T                   = cluster(dendrogram,'maxclust',N_clusters_selected);
                    
                    hh = figure();
                    hold on
                    plot([1 N_c],DB_vec,'Color','b','Linestyle','-','Linewidth',2,'Marker','s','Markersize',6)
                    plot([1 N_c],DB_av_vec,'Color','r','Linestyle','-','Linewidth',2,'Marker','*','Markersize',6)
                    plot([1 N_c],s_mean_vec,'Color','k','Linestyle','-','Linewidth',2,'Marker','o','Markersize',6)
                    plot([1 N_c],s_av_mean_vec,'Color','m','Linestyle','-','Linewidth',2,'Marker','>','Markersize',6)
                    h_legend = legend('Davies-Bouldin Index','Weighted Davies-Bouldin Index',...
                        'Average Silhouette','Weighted Average Silhouette','Location','Best');
                    set(h_legend,'Fontname','Avantgarde','Fontsize',10)
                    xlabel('Number of clusters','Fontname','Avantgarde','Fontsize',14)
                    ylabel('Figure of merit','Fontname','Avantgarde','Fontsize',14)
                    title([num2str(airport_ID(i)),'-',num2str(other_airports(j)),...
                          ': ',num2str(N_traj),' trajectories ',num2str(N_clusters_selected),' clusters - ',...
                          'dist airp ',num2str(dist_airp),' - ',' av. FD ',num2str(mean_FD)],'Fontname','Avantgarde','Fontsize',14)
                    grid on
                    
                    saveas(hh,[strcat(pwd,'/FIGURES/',num2str(airport_ID(i)),'_',num2str(other_airports(j)),'.eps')],'epsc')
                    
                end
                
                % Building matrix with flights belonging to the same
                % cluster along each row
                CLUSTERS = zeros(N_clusters_selected,N_traj);
                for ii=1:N_clusters_selected
                    idx = find(T==ii);
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
            end
            
            disp(['Number of clusters: ',num2str(N_clusters_selected)]);
            
        else
        end
        disp('%%%%%%%%%%%%%%%%%%%%%');
    end
end