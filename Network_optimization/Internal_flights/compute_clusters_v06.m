clc
clear all
close all

% Checking if the folder exists already or not
if isequal(exist(strcat(pwd,'/NETWORK/INT_Clusters'),'dir'),7) % 7 = directory
else
    mkdir(strcat(pwd,'/NETWORK/INT_Clusters'));
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
            c_lim     = 0.3;
            
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
            %disp(['Number of outliers removed: ',num2str(numel(outliers))]);
            %disp(['Number of trajectories remaining: ',num2str(numel(flights))]);
            
            %dlmwrite((strcat(pwd,'/NETWORK/INT_FDMatrix/',num2str(airport_ID(i)),...
            %    '/',num2str(airport_ID(i)),'_',num2str(other_airports(j)),'_red.txt')),FM,' ');
            %dlmwrite((strcat(pwd,'/NETWORK/INT_FDMatrix/',num2str(airport_ID(i)),...
            %    '/',num2str(airport_ID(i)),'_',num2str(other_airports(j)),'_flightsID.txt')),flights,' ');
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Selection of the optimal number of clusters %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            N_traj = numel(flights);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            dist_mat_km              = FM+FM';
            dist_vec_km              = squareform(dist_mat_km);
            tree                     = linkage(dist_vec_km,'average');
            threshold_single_cluster = 30; % [km]
            
            % If this condition is satisfied, it means we have a very
            % compact set of trajectories. Thus, the best solution is to
            % consider a single cluster
            if tree(end,3)<=threshold_single_cluster
                NC_opt   = 1;
                CLUSTERS = flights';
                disp(['Number of clusters from airport ',num2str(airport_ID(i)),' to airport ',num2str(other_airports(j)),' is 1: very compact cluster'])
            else
                LB          = 50;  % [km]
                UB          = 200; % [km]
                idx_LB      = find(tree(:,3)>=LB);
                idx_UB      = find(tree(:,3)<=UB);
                idx         = intersect(idx_LB,idx_UB);
                
                if isempty(idx)
                    idx = [N_traj-3;N_traj-2;N_traj-1];
                else
                end
                
                N_c         = sort(N_traj-idx);
                max_N_c     = 20;
                
                if N_c(1) == 1
                    N_c = [N_c(2:end);N_c(end)+1];
                    % At most, N_c(end) = N_traj-1, thus N_c(end)+1
                    % will be N_traj in the worst case scenario (it
                    % corresponds to the case where each singleton
                    % forms its own cluster)
                else
                    N_c = [N_c;N_c(end)+1];
                end
                
                if numel(N_c)>max_N_c
                    N_c = N_c(1:max_N_c);
                else
                end
                
                if numel(N_c)==1
                    if N_c == 2
                        N_c = [2;3];
                    elseif N_c == N_traj-1
                        N_c = [N_traj-2;N_traj-1];
                    else
                        N_c = [N_c;N_c+1];
                    end
                else
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Using hierarchical clustering to compute the number of clusters %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                AS   = zeros(numel(N_c),1);
                DB   = zeros(numel(N_c),1);
                DI   = zeros(numel(N_c),1);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Computing the 3 performance indices %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                tree         = linkage(dist_vec_km,'average');
                flag_indices = 0;
                
                for k=1:numel(N_c)
                    
                    N_clust    = N_c(k);
                    T          = cluster(tree,'maxclust',N_clust);
                    
                    s      = silhouette([],T,dist_vec_km);
                    AS(k)  = mean(s);
                    DB(k)  = DaviesBouldin(dist_mat_km,N_clust,T,flag_indices);
                    DI(k)  = Dunn(dist_mat_km,N_clust,T,flag_indices);
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Extending the N_c vector with an additional element before  %%%
                %%% N_c(1), for the same reason already introduced for          %%%
                %%% N_c(end). If N_c(1)=2, then we artificially add values for  %%%
                %%% the three indices (they are not defined for N_c=1),         %%%
                %%% otherwise we compute the values                             %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
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
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Computing the optimal number of clusters comparing %%%
                %%% the performances of the different indices          %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % We use findpeaks as provided for Average Silhouette and
                % Dunn Index (we want to maximize these indices), while we
                % add a minus in front of the Davies-Bouldin Index (because
                % we want to minimize it)
                
                local_maxmin    = zeros(3,1);
                NC_optimal      = NaN(3,1);
                
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
                
                Nc_index1  = NC_optimal(1);
                Nc_index2  = NC_optimal(2);
                Nc_index3  = NC_optimal(3);
                NC_opt = min([Nc_index1,Nc_index2,Nc_index3]);
                
                if isnan(NC_opt)
                    NC_opt = 2;
                else
                end
                
                disp(['Number of clusters from airport ',num2str(airport_ID(i)),' to airport ',num2str(other_airports(j)),' is: ',num2str(NC_opt)])
                
            end
            
            % Storing matrix in the associated folder
            %dlmwrite(strcat(pwd,'/NETWORK/INT_Clusters/',num2str(airport_ID(i)),'/',num2str(airport_ID(i)),'_',num2str(other_airports(j)),'.txt'),CLUSTERS,' ');
            %disp(['Number of clusters: ',num2str(NC_opt)]);
            
        else
            disp(['No valid trajectories from ',num2str(airport_ID(i)),' to airport ',num2str(other_airports(j))])
        end
    end
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
end

            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
%             
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             
%             % If we have only 2 trajectories characterizing this OD pair,
%             % we check whether the onlt non-zero component of the Frechet
%             % Distance Matrix is bigger or smaller than a specific
%             % threshold. In the first case, we consider each trajectory as
%             % a single cluster. In the second case, the two trajectories
%             % form a unique cluster
%             if N_traj<5
%                 NC_opt   = 1;
%                 CLUSTERS = flights';
%             else
%                 
%                 
%                 % Reshape dist_mat_km into a vector compatible with the matlab's pdist
%                 % format
%                 dist_vec_km = [];
%                 for ii=1:N_traj-1
%                     dist_vec_km= horzcat(dist_vec_km,FM(ii,ii+1:end));
%                 end
%                 
%                 max_clusters = 8;
%                 N_c          = 2:min(N_traj-1,max_clusters);
%                 
%                 R_e         = 6371;  % Earth's radius [km]
%                 h           = 10;    % altitude [km]
%                 radius      = R_e+h; % radius of the sphere where we project the points
%                 dist_mat_km = FM;
%                 
%                 AS   = zeros(numel(N_c),1);
%                 DB   = zeros(numel(N_c),1);
%                 DI   = zeros(numel(N_c),1);
%                 
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 %%% Computing the 3 performance indices %%%
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 
%                 for k=1:numel(N_c)
%                     
%                     N_clust    = N_c(k);
%                     dendrogram = linkage(dist_vec_km,'average');
%                     T          = cluster(dendrogram,'maxclust',N_clust);
%                     
%                     s      = silhouette([],T,dist_vec_km);
%                     AS(k)  = mean(s);
%                     DB(k)  = DaviesBouldin(FM,N_clust,T);
%                     DI(k)  = Dunn(FM,N_clust,T);
%                     
%                 end
%                 
%                 N_c_aug = [1 N_c];
%                 AS_aug  = [0;AS];
%                 DB_aug  = [1;DB];
%                 DI_aug  = [1;DI];
%                 
%                 
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 %%% Computing the optimal number of clusters comparing %%%
%                 %%% the performances of the different indices          %%%
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 
%                 % We use findpeaks as provided for Average Silhouette and
%                 % Dunn Index (we want to maximize these indices), while we
%                 % add a minus in front of the Davies-Bouldin Index (because
%                 % we want to minimize it)
%                 
%                 local_maxmin    = zeros(3,1);
%                 NC_optimal      = zeros(3,1);
%                 
%                 [val_AS,idx_AS] = findpeaks(AS_aug);
%                 if isempty(val_AS)
%                 else
%                     local_maxmin(1,1) = 1;
%                     if numel(val_AS(:,1))>1
%                         [~,dummy]       = max(val_AS(:,1));
%                         NC_optimal(1,1) = N_c_aug(idx_AS(dummy));
%                     else
%                         NC_optimal(1,1) = N_c_aug(idx_AS);
%                     end
%                 end
%                 [val_DB,idx_DB] = findpeaks(-DB_aug);
%                 if isempty(val_DB)
%                 else
%                     local_maxmin(2,1) = 1;
%                     if numel(val_DB(:,1))>1
%                         [~,dummy]       = max(val_DB(:,1));
%                         NC_optimal(2,1) = N_c_aug(idx_DB(dummy));
%                     else
%                         NC_optimal(2,1) = N_c_aug(idx_DB);
%                     end
%                 end
%                 [val_DI,idx_DI] = findpeaks(DI_aug);
%                 if isempty(val_DI)
%                 else
%                     local_maxmin(3,1) = 1;
%                     if numel(val_DI(:,1))>1
%                         [~,dummy]       = max(val_DI(:,1));
%                         NC_optimal(3,1) = N_c_aug(idx_DI(dummy));
%                     else
%                         NC_optimal(3,1) = N_c_aug(idx_DI);
%                     end
%                 end
%                 
%                 % if none of the three indices displays a local max/min,
%                 % then we assume the trajectories form a single cluster
%                 if isempty(find(local_maxmin==1,1))
%                     NC_opt = 1;
%                 else
%                     all_values     = [AS_aug DB_aug DI_aug];
%                     N_good_indices = numel(find(local_maxmin==1));
%                     % only one index displays a local max/min. The
%                     % associated number of clusters if the optimal number
%                     % of clusters
%                     if N_good_indices == 1
%                         good_index = find(local_maxmin==1);
%                         NC_opt     = NC_optimal(good_index);
%                         % two indices display a local max/min. We pick the
%                         % number of clusters that introduces the smallest loss
%                         % of performance in the other index
%                     elseif N_good_indices == 2
%                         good_index = find(local_maxmin==1);
%                         index1     = good_index(1);
%                         Nc_index1  = NC_optimal(index1);
%                         index2     = good_index(2);
%                         Nc_index2  = NC_optimal(index2);
%                         % evaluating loss of performance of first index
%                         % when using optimal number of clusters of the
%                         % second index
%                         lop1       = abs(abs(all_values(Nc_index1,index1)-all_values(Nc_index2,index1))/all_values(Nc_index1,index1)*100);
%                         % evaluating loss of performance of second index
%                         % when using optimal number of clusters of the
%                         % first index
%                         lop2       = abs(abs(all_values(Nc_index2,index2)-all_values(Nc_index1,index2))/all_values(Nc_index2,index2)*100);
%                         
%                         if lop1 < lop2
%                             NC_opt = Nc_index2;
%                         else
%                             NC_opt = Nc_index1;
%                         end
%                         % all three indices display a local max/min. We pick the
%                         % number of clusters that introduces the smallest loss
%                         % of performance when applied to the other two indices
%                     else
%                         Nc_index1  = NC_optimal(1);
%                         Nc_index2  = NC_optimal(2);
%                         Nc_index3  = NC_optimal(3);
%                         % evaluating average loss of performance of
%                         % second/third indice when using optimal number of
%                         % clusters of the first index
%                         lop1       = (abs(abs(all_values(Nc_index1,2)-all_values(Nc_index2,2))/all_values(Nc_index2,2)*100)+...
%                             abs(abs(all_values(Nc_index1,3)-all_values(Nc_index3,3))/all_values(Nc_index3,3)*100))/2;
%                         % evaluating average loss of performance of
%                         % first/third indice when using optimal number of
%                         % clusters of the second index
%                         lop2       = (abs(abs(all_values(Nc_index2,1)-all_values(Nc_index1,1))/all_values(Nc_index1,1)*100)+...
%                             abs(abs(all_values(Nc_index2,3)-all_values(Nc_index3,3))/all_values(Nc_index3,3)*100))/2;
%                         % evaluating average loss of performance of
%                         % first/second indice when using optimal number of
%                         % clusters of the third index
%                         lop3       = (abs(abs(all_values(Nc_index3,1)-all_values(Nc_index1,1))/all_values(Nc_index1,1)*100)+...
%                             abs(abs(all_values(Nc_index3,2)-all_values(Nc_index2,2))/all_values(Nc_index2,2)*100))/2;
%                         
%                         Nc_all  = [Nc_index1 Nc_index2 Nc_index3];
%                         [~,idx] = min([lop1 lop2 lop3]);
%                         NC_opt  = Nc_all(idx);
%                         
%                     end
%                 end
%                 
%                 % Building matrix with flights belonging to the same
%                 % cluster along each row
%                 T        = cluster(dendrogram,'maxclust',NC_opt);
%                 CLUSTERS = zeros(NC_opt,N_traj);
%                 for ii=1:NC_opt
%                     idx = find(T==ii);
%                     if ii==1
%                         max_traj = numel(idx);
%                     else
%                         if numel(idx)>max_traj
%                             max_traj = numel(idx);
%                         else
%                         end
%                     end
%                     CLUSTERS(ii,1:numel(idx)) = transpose(flights(idx));
%                 end
%                 % Getting rid of extra columns
%                 CLUSTERS = CLUSTERS(:,1:max_traj);
%                 
%             end
            
            

