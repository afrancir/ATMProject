% Testing performances of different indices when assessing the optimal
% number of clusters

clc
clear all
close all

Origin      = 33;
Destination = 76;

% Loading all trajectories
ALL_TRAJ = load(strcat(pwd,'/NETWORK/INT_traj_filt_info/',num2str(Origin),'/',num2str(Origin),'_',num2str(Destination),'_int.txt'));
flights = unique(ALL_TRAJ(:,1)); % number of trajectories

% Loading Frechet Distance Matrix
FM = load(strcat(pwd,'/NETWORK/INT_FDMatrix/',num2str(Origin),...
    '/',num2str(Origin),'_',num2str(Destination),'.txt'));

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Selection of the optimal number of clusters %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_traj = numel(flights);

% If we have only 2 trajectories characterizing this OD pair,
% we check whether the onlt non-zero component of the Frechet
% Distance Matrix is bigger or smaller than a specific
% threshold. In the first case, we consider each trajectory as
% a single cluster. In the second case, the two trajectories
% form a unique cluster
if N_traj==2
    origin_airp_idx = find(AIRPORTS(:,1)==Origin);
    dest_airp_idx   = find(AIRPORTS(:,1)==Destination);
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
    
    
    % Reshape dist_mat_km into a vector compatible with the matlab's pdist
    % format
    dist_vec_km = [];
    for ii=1:N_traj-1
        dist_vec_km= horzcat(dist_vec_km,FM(ii,ii+1:end));
    end
    
    max_clusters = 8;
    N_c          = 2:min(N_traj,max_clusters);
    
    R_e         = 6371;  % Earth's radius [km]
    h           = 10;    % altitude [km]
    radius      = R_e+h; % radius of the sphere where we project the points
    dist_mat_km = FM;
    
    FM_full     = FM+FM';
    
    ASW   = zeros(numel(N_c),1);
    DUNN  = zeros(numel(N_c),1);
    DB    = zeros(numel(N_c),1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Computing the 4 performance indexes %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for j=1:numel(N_c)
        Nc = N_c(j);
        dendrogram = linkage(dist_vec_km,'average');
        T          = cluster(dendrogram,'maxclust',Nc);
        s          = silhouette([],T,dist_vec_km);
        ASW(j)     = mean(s);
        DUNN(j)    = Dunn(FM,Nc,T);
        DB(j)      = DaviesBouldin(FM,Nc,T);
    end
    
    [~,idx_ASW]  = max(ASW);
    [~,idx_DUNN] = max(DUNN);
    [~,idx_DB]   = min(DB);
    
    N_c(idx_ASW)
    N_c(idx_DUNN)
    N_c(idx_DB)
    
    
    figure()
    hold on
    plot(N_c,ASW,'Color','b','Linestyle','-','Linewidth',2,'Marker','s','Markersize',6)
    plot(N_c,DUNN,'Color','r','Linestyle','-','Linewidth',2,'Marker','*','Markersize',6)
    plot(N_c,DB,'Color','k','Linestyle','-','Linewidth',2,'Marker','o','Markersize',6)
    h_legend = legend('Average Silhouette','Dunn Index','Davies-Bouldin Index'...
        ,'Location','Best');
    set(h_legend,'Fontname','Avantgarde','Fontsize',10)
    xlabel('Number of clusters','Fontname','Avantgarde','Fontsize',14)
    ylabel('Figure of merit','Fontname','Avantgarde','Fontsize',14)
    grid on
    
end
             